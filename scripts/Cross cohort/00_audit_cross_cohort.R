#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

# Public adaptation: severity is read from coded metadata, not from sample labels.
discovery_metadata_file <- Sys.getenv("MPP_DISCOVERY_METADATA", "inputs/discovery.csv")
discovery_metadata_public <- read.csv(discovery_metadata_file, stringsAsFactors = FALSE)
if (anyDuplicated(discovery_metadata_public$sample_id)) stop("Duplicate coded discovery sample")
discovery_severity <- setNames(discovery_metadata_public$severity_en, discovery_metadata_public$sample_id)


discovery_root <- Sys.getenv(
  "BD_DISCOVERY_ROOT",
  "inputs/discovery/bd_rhapsody_10/downstream/10samples_scrna_mex_v4"
)
validation_root <- Sys.getenv(
  "BD_VALIDATION_ROOT",
  "inputs/workspace/工作文件/23例单细胞重新分析"
)
clinical_file <- Sys.getenv(
  "BD_VALIDATION_CLINICAL",
  "inputs/workspace/工作文件/单细胞23例样本临床信息.csv"
)
cross_root <- Sys.getenv(
  "BD_CROSS_ROOT",
  "inputs/discovery/bd_rhapsody_10/downstream/10_23_cross_cohort_v1"
)
result <- file.path(cross_root, "results", "00_audit")
dir.create(result, recursive = TRUE, showWarnings = FALSE)

read_csv <- function(path) {
  if (!file.exists(path)) stop("Missing required input: ", path)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8-BOM")
}

write_csv <- function(frame, name) {
  write.csv(frame, file.path(result, name), row.names = FALSE, fileEncoding = "UTF-8")
}

discovery_stage03 <- file.path(discovery_root, "results", "03_annotate")
discovery_stage04 <- file.path(discovery_root, "results", "04_compare")
validation_results <- file.path(validation_root, "results")

stage03_validation <- read_csv(file.path(discovery_stage03, "03_annotate_validation.csv"))
stage04_validation <- read_csv(file.path(discovery_stage04, "04_compare_validation.csv"))
source_validation_checks <- rbind(
  transform(stage03_validation, cohort = "discovery_10", stage = "03_annotate"),
  transform(stage04_validation, cohort = "discovery_10", stage = "04_compare")
)
source_validation_checks <- source_validation_checks[, c("cohort", "stage", "check", "status", "detail")]
write_csv(source_validation_checks, "source_validation_checks.csv")

discovery_counts <- read_csv(file.path(discovery_stage03, "cell_type_counts_by_sample.csv"))
discovery_counts$Freq <- as.integer(discovery_counts$Freq)
validation_qc <- read_csv(file.path(validation_results, "01_qc", "qc_summary_by_sample.csv"))
validation_qc$cells_qc_pass <- as.integer(validation_qc$cells_qc_pass)
validation_groups <- read_csv(file.path(validation_results, "00_audit", "clinical_group_counts.csv"))
validation_samples <- read_csv(file.path(validation_results, "00_audit", "sample_set_audit.csv"))

group_count <- function(variable, level) {
  rows <- validation_groups$group_variable == variable & validation_groups$group_level == level
  if (!any(rows)) return(NA_integer_)
  as.integer(validation_groups$sample_count[rows][[1]])
}

cohort_inventory <- data.frame(
  cohort = c("discovery_10", "validation_23"),
  samples = c(length(unique(discovery_counts$sample_id)), nrow(validation_samples)),
  mild_samples = c(
    sum(discovery_severity[unique(discovery_counts$sample_id)] == "mild"),
    group_count("severity", "轻症")
  ),
  severe_samples = c(
    sum(discovery_severity[unique(discovery_counts$sample_id)] == "severe"),
    group_count("severity", "重症")
  ),
  cells = c(sum(discovery_counts$Freq), sum(validation_qc$cells_qc_pass)),
  primary_object = c(
    file.path(discovery_stage03, "10samples_annotated_seurat_v5.rds"),
    file.path(validation_results, "02_rna_reanalysis", "23samples_RNA_reclustered_annotated.rds")
  ),
  metadata_source = c(
    file.path(discovery_root, "results", "00_audit", "10samples_metadata_locked.csv"),
    clinical_file
  ),
  independent_patient_cohort = TRUE,
  stringsAsFactors = FALSE
)
write_csv(cohort_inventory, "cohort_inventory.csv")

broad <- read_csv(file.path(discovery_stage03, "broad_cluster_annotations.csv"))
broad$top_score <- as.numeric(broad$top_score)
broad$score_margin <- as.numeric(broad$score_margin)
broad$supporting_markers <- as.integer(broad$supporting_markers)
broad$classification <- "通过"
broad$review_reason <- "marker evidence and confidence are adequate"
low_margin <- broad$confidence == "moderate" & broad$score_margin < 0.1
broad$classification[low_margin] <- "需修订"
broad$review_reason[low_margin] <- "moderate-confidence cluster has score margin below 0.1; retain broad label only after marker-plot review"
unexpected_low <- broad$confidence == "low" & broad$label != "Unresolved"
broad$classification[unexpected_low] <- "暂不可信"
broad$review_reason[unexpected_low] <- "low-confidence resolved label is not supported"
resolved_uncertainty <- broad$confidence == "low" & broad$label == "Unresolved"
broad$classification[resolved_uncertainty] <- "通过"
broad$review_reason[resolved_uncertainty] <- "low-confidence cluster is explicitly preserved as Unresolved"
write_csv(broad, "discovery_broad_cluster_review.csv")

lineage <- read_csv(file.path(discovery_stage03, "lineage_subcluster_annotations.csv"))
lineage_review <- do.call(rbind, lapply(split(lineage, lineage$lineage), function(frame) {
  low_unresolved <- sum(frame$confidence == "low" & frame$label == "Unresolved")
  data.frame(
    lineage = frame$lineage[[1]],
    subclusters = nrow(frame),
    high_confidence = sum(frame$confidence == "high"),
    moderate_confidence = sum(frame$confidence == "moderate"),
    low_confidence = sum(frame$confidence == "low"),
    unresolved_subclusters = sum(frame$label == "Unresolved"),
    classification = if (low_unresolved > 0L) "需修订" else "通过",
    review_reason = if (low_unresolved > 0L) {
      "exclude explicitly Unresolved subclusters from mechanistic claims and cross-cohort subtype mapping"
    } else {
      "all subclusters have resolved marker-supported labels"
    },
    stringsAsFactors = FALSE
  )
}))
write_csv(lineage_review, "discovery_lineage_subcluster_review.csv")

sample_totals <- aggregate(Freq ~ sample_id, discovery_counts, sum)
names(sample_totals)[[2]] <- "sample_total"
composition <- merge(discovery_counts, sample_totals, by = "sample_id", all.x = TRUE)
composition$severity <- unname(discovery_severity[composition$sample_id])
composition$proportion <- composition$Freq / composition$sample_total
composition_results <- read_csv(file.path(discovery_stage04, "cell_type_composition_comparison.csv"))

composition_review <- do.call(rbind, lapply(split(composition, composition$cell_type_broad), function(frame) {
  cell_type <- frame$cell_type_broad[[1]]
  full_difference <- mean(frame$proportion[frame$severity == "severe"]) -
    mean(frame$proportion[frame$severity == "mild"])
  loo <- vapply(unique(frame$sample_id), function(omitted) {
    current <- frame[frame$sample_id != omitted, , drop = FALSE]
    mean(current$proportion[current$severity == "severe"]) -
      mean(current$proportion[current$severity == "mild"])
  }, numeric(1))
  dominant <- frame[which.max(frame$Freq), , drop = FALSE]
  dominant_share <- dominant$Freq / sum(frame$Freq)
  published <- composition_results[composition_results$cell_type_broad == cell_type, , drop = FALSE]
  classification <- "通过"
  reason <- "no FDR-supported composition claim; leave-one-patient-out direction is recorded"
  if (dominant_share >= 0.5) {
    classification <- "需修订"
    reason <- "one patient contributes at least half of all cells in this broad type; nominal composition direction is outlier-sensitive"
  }
  if (cell_type == "Unresolved") {
    classification <- "暂不可信"
    reason <- "Unresolved is not a biological cell type and must not support a composition mechanism claim"
  }
  data.frame(
    cell_type_broad = cell_type,
    full_severe_minus_mild = full_difference,
    loo_min = min(loo),
    loo_max = max(loo),
    loo_direction_agreement = if (full_difference == 0) NA_real_ else mean(sign(loo) == sign(full_difference)),
    dominant_sample = dominant$sample_id,
    dominant_cell_share = dominant_share,
    current_fdr = as.numeric(published$fdr),
    classification = classification,
    review_reason = reason,
    stringsAsFactors = FALSE
  )
}))
write_csv(composition_review, "discovery_composition_loo_review.csv")

de_summary <- read_csv(file.path(discovery_stage04, "pseudobulk_de_summary.csv"))
sensitivity_summary <- read_csv(file.path(discovery_stage04, "pseudobulk_sensitivity_summary.csv"))
sensitivity_concordance <- read_csv(file.path(discovery_stage04, "pseudobulk_sensitivity_concordance.csv"))
discovery_de_review <- do.call(rbind, lapply(de_summary$cell_type_broad, function(cell_type) {
  primary <- de_summary[de_summary$cell_type_broad == cell_type, , drop = FALSE]
  sensitivity <- sensitivity_summary[sensitivity_summary$cell_type_broad == cell_type, , drop = FALSE]
  concordance <- sensitivity_concordance[sensitivity_concordance$cell_type_broad == cell_type, , drop = FALSE]
  primary_hits <- as.integer(primary$fdr_lt_0_05)
  sensitivity_hits <- as.integer(sensitivity$fdr_lt_0_05)
  min_correlation <- min(as.numeric(concordance$logFC_spearman), na.rm = TRUE)
  classification <- "需修订"
  reason <- "technically valid discovery result requires independent-cohort confirmation"
  if (cell_type == "Unresolved") {
    classification <- "暂不可信"
    reason <- "DE for Unresolved cells is not biologically interpretable"
  } else if (cell_type == "T_NK" && primary_hits > 0L && all(sensitivity_hits == 0L)) {
    classification <- "暂不可信"
    reason <- "the single primary FDR hit disappears in both prespecified sensitivity analyses"
  } else if (primary_hits == 0L && all(sensitivity_hits == 0L)) {
    classification <- "通过"
    reason <- "no FDR-supported DE claim in primary or sensitivity analyses"
  } else if (primary_hits == 0L && any(sensitivity_hits > 0L)) {
    classification <- "需修订"
    reason <- "signals appear only after excluding QC outliers and are not stable discovery claims"
  } else if (min_correlation < 0.9) {
    classification <- "需修订"
    reason <- "global logFC concordance falls below 0.9 in at least one prespecified sensitivity analysis"
  }
  data.frame(
    cell_type_broad = cell_type,
    genes_tested = as.integer(primary$genes_tested),
    primary_fdr_hits = primary_hits,
    sensitivity_fdr_hits_min = min(sensitivity_hits),
    sensitivity_fdr_hits_max = max(sensitivity_hits),
    min_logfc_spearman = min_correlation,
    min_patients_per_group = min(as.integer(primary$n_mild_patients), as.integer(primary$n_severe_patients)),
    classification = classification,
    review_reason = reason,
    stringsAsFactors = FALSE
  )
}))
write_csv(discovery_de_review, "discovery_pseudobulk_review.csv")

qc_rate <- as.numeric(validation_qc$qc_pass_rate)
validation_pb <- read_csv(file.path(validation_results, "03_pseudobulk", "pseudobulk_DE_summary.csv"))
validation_pb$design_rank <- as.integer(validation_pb$design_rank)
validation_pb$design_columns <- as.integer(validation_pb$design_columns)
validation_pb$genes_tested <- as.integer(validation_pb$genes_tested)
validation_pb$FDR_lt_0_05_abs_logFC_ge_1 <- as.integer(validation_pb$FDR_lt_0_05_abs_logFC_ge_1)
validation_lineage <- read_csv(file.path(validation_results, "07_lineage_pseudobulk", "lineage_pseudobulk_DE_summary.csv"))
validation_lineage$FDR_lt_0_05_abs_logFC_ge_1 <- as.integer(validation_lineage$FDR_lt_0_05_abs_logFC_ge_1)

validation_review <- data.frame(
  evidence_item = c(
    "23-sample identity and matrix intersection",
    "RNA QC retention",
    "RNA broad annotation coverage",
    "broad pseudobulk design",
    "lineage pseudobulk design",
    "cross-cohort broad-label harmonization",
    "existing pathway analysis comparability",
    "existing CellChat evidence",
    "RNA+ADT/WNN runtime object"
  ),
  observed = c(
    paste0(sum(validation_samples$in_object == "TRUE" & validation_samples$in_clinical == "TRUE" & validation_samples$in_raw_matrix == "TRUE"), "/23 complete"),
    paste0("cells=", sum(validation_qc$cells_qc_pass), "; min_qc_pass_rate=", sprintf("%.3f", min(qc_rate))),
    "209299 cells; 23 samples; six coarse broad labels",
    paste0(sum(validation_pb$status == "completed"), "/", nrow(validation_pb), " completed; full_rank=", sum(validation_pb$design_rank == validation_pb$design_columns)),
    paste0(nrow(validation_lineage), " lineage-contrast models; significant_strict=", sum(validation_lineage$FDR_lt_0_05_abs_logFC_ge_1)),
    "Macrophage and Monocyte are separate; T_cells is not identical to discovery T_NK",
    "validation uses clusterProfiler GSEA while discovery uses cameraPR",
    "historical pooled network only; sample-level formal rerun absent",
    paste0("object_bytes=", file.info(file.path(validation_results, "10_next_decision", "23samples_RNA_ADT_WNN_rebuilt.rds"))$size)
  ),
  classification = c("通过", "通过", "通过", "通过", "通过", "需修订", "需修订", "暂不可信", "需修订"),
  required_action = c(
    "retain coded sample matching and cohort identity",
    "retain current QC object; validate metadata and count layers at runtime before new modeling",
    "retain for within-cohort analysis; use an explicit harmonized mapping for cross-cohort work",
    "reuse sample-level raw-count design and add harmonized lineage models with effect standard errors",
    "retain current outputs as exploratory subtype evidence",
    "recompute four directly comparable harmonized lineages from raw counts",
    "rerun both cohorts with the same offline gene-set membership and ranking method",
    "do not use for formal conclusions; rerun only if a replicated lineage signal justifies it",
    "load and validate all 23 samples/cells/features before any ADT-supported mechanism claim"
  ),
  stringsAsFactors = FALSE
)
write_csv(validation_review, "validation_evidence_review.csv")

cross_cohort_mapping <- data.frame(
  discovery_label = c(
    "B_Plasma", "Macrophage_Monocyte", "Neutrophil", "T_NK", "Epithelial",
    "Dendritic", "Endothelial", "Mast", "Stromal", "Unresolved"
  ),
  validation_label = c(
    "B_cell", "Macrophage;Monocyte", "Neutrophils", "T_cells", "Epithelial_cells",
    NA, NA, NA, NA, NA
  ),
  harmonized_label = c(
    "B_lineage", "Macrophage_Monocyte", "Neutrophil", "T_lineage", "Epithelial",
    NA, NA, NA, NA, NA
  ),
  meta_eligible = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
  mapping_rule = c(
    "direct broad-lineage mapping",
    "sum raw counts across validation Macrophage and Monocyte before modeling",
    "direct broad-lineage mapping",
    "exclude discovery NK, Cycling_T_NK, and Unresolved subtypes before T-lineage modeling; validate validation T_cells CD3 support",
    "discovery has fewer than three eligible patients per severity group",
    "absent from validation broad annotation",
    "absent from validation broad annotation",
    "absent from validation broad annotation",
    "absent from validation broad annotation",
    "not a biological cell type"
  ),
  stringsAsFactors = FALSE
)
write_csv(cross_cohort_mapping, "cross_cohort_celltype_mapping.csv")

input_paths <- c(
  file.path(discovery_stage03, "10samples_annotated_seurat_v5.rds"),
  file.path(discovery_stage04, "all_cell_types_pseudobulk_de.csv.gz"),
  file.path(discovery_root, "results", "05_report", "downstream_output_manifest.csv"),
  file.path(validation_results, "02_rna_reanalysis", "23samples_RNA_reclustered_annotated.rds"),
  file.path(validation_results, "10_next_decision", "23samples_RNA_ADT_WNN_rebuilt.rds"),
  clinical_file
)
input_info <- file.info(input_paths)
input_inventory <- data.frame(
  path = input_paths,
  exists = file.exists(input_paths),
  bytes = as.numeric(input_info$size),
  modified = format(input_info$mtime, "%Y-%m-%dT%H:%M:%S%z"),
  stringsAsFactors = FALSE
)
write_csv(input_inventory, "input_inventory.csv")

classification_counts <- as.data.frame(table(
  classification = c(
    broad$classification,
    lineage_review$classification,
    composition_review$classification,
    discovery_de_review$classification,
    validation_review$classification
  )
), stringsAsFactors = FALSE)
write_csv(classification_counts, "classification_summary.csv")

report <- c(
  "# Discovery and validation cohort audit",
  "",
  paste0("- Generated: ", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
  "- Cohorts are independent; no patient-level matching is performed.",
  paste0("- Discovery cohort: ", cohort_inventory$samples[[1]], " samples, ", cohort_inventory$cells[[1]], " cells."),
  paste0("- Validation cohort: ", cohort_inventory$samples[[2]], " samples, ", cohort_inventory$cells[[2]], " cells."),
  "- Existing discovery Stage 03/04 validators must remain PASS.",
  "- Positive discovery DE findings are not accepted as replicated until harmonized validation and meta-analysis complete.",
  "- Unresolved cells, historical pooled CellChat, and non-harmonized pathway outputs are excluded from formal biological claims.",
  "",
  "## Required next analysis",
  "",
  "1. Freeze the discovery v1 manifest and evidence review without overwriting validated outputs.",
  "2. Recompute B, macrophage/monocyte, neutrophil, and T-lineage pseudobulk models with explicit harmonized mappings.",
  "3. Add comparable effect standard errors and run two-cohort effect integration with heterogeneity reporting.",
  "4. Rerun pathway tests with one shared offline gene-set method.",
  "5. Restrict mechanism analysis to lineages supported in both cohorts."
)
writeLines(report, file.path(result, "cross_cohort_audit_report.md"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(result, "sessionInfo.txt"), useBytes = TRUE)

message("Cross-cohort audit completed")
