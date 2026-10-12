#!/usr/bin/env Rscript

lib <- Sys.getenv("BD_R_LIB", "inputs/discovery/bd_rhapsody_10/downstream/10samples_scrna/env/host-Rlib")
.libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(edgeR)
  library(limma)
})

options(stringsAsFactors = FALSE, warn = 1)
set.seed(20260727)

cross_root <- Sys.getenv(
  "BD_CROSS_ROOT",
  "inputs/discovery/bd_rhapsody_10/downstream/10_23_cross_cohort_v1"
)
audit_root <- file.path(cross_root, "results", "00_audit")
freeze_root <- file.path(cross_root, "results", "01_discovery_freeze")
result <- file.path(cross_root, "results", "02_harmonized_validation")
count_root <- file.path(result, "pseudobulk_counts")
dir.create(count_root, recursive = TRUE, showWarnings = FALSE)

read_csv_path <- function(path) {
  if (!file.exists(path)) stop("Missing required input: ", path)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8-BOM")
}

write_csv <- function(frame, name) {
  write.csv(frame, file.path(result, name), row.names = FALSE, fileEncoding = "UTF-8")
}

freeze <- read_csv_path(file.path(freeze_root, "discovery_v1_freeze_record.csv"))
if (nrow(freeze) != 1L || freeze$freeze_id != "discovery_10_v1") stop("Discovery v1 freeze record is invalid")
inventory <- read_csv_path(file.path(audit_root, "cohort_inventory.csv"))
discovery_object_path <- inventory$primary_object[inventory$cohort == "discovery_10"]
validation_object_path <- inventory$primary_object[inventory$cohort == "validation_23"]
if (length(discovery_object_path) != 1L || length(validation_object_path) != 1L) stop("Cohort inventory is incomplete")
validation_root <- sub("/results/.*$", "", validation_object_path)
validation_count_path <- file.path(validation_root, "results", "03_pseudobulk", "pseudobulk_raw_counts.rds")
validation_meta_path <- file.path(validation_root, "results", "03_pseudobulk", "pseudobulk_sample_metadata.csv")

lineages <- c("B_lineage", "Macrophage_Monocyte", "Neutrophil", "T_lineage")
discovery_broad_map <- c(
  B_Plasma = "B_lineage",
  Macrophage_Monocyte = "Macrophage_Monocyte",
  Neutrophil = "Neutrophil",
  T_NK = "T_lineage"
)
validation_label_map <- list(
  B_lineage = "B_cell",
  Macrophage_Monocyte = c("Macrophage", "Monocyte"),
  Neutrophil = "Neutrophils",
  T_lineage = "T_cells"
)
forbidden_t_subtypes <- c("NK", "Cycling_T_NK", "Unresolved")

aggregate_cells <- function(counts, cell_meta, group_label, sample_info, total_cells, cell_weights = NULL) {
  if (is.null(cell_weights)) cell_weights <- rep(1, nrow(cell_meta))
  if (length(cell_weights) != nrow(cell_meta)) stop("Cell-count weights do not match aggregation metadata")
  keep <- !is.na(group_label)
  selected_meta <- cell_meta[keep, , drop = FALSE]
  selected_meta$harmonized_lineage <- group_label[keep]
  selected_meta$group_id <- paste(selected_meta$sample_id, selected_meta$harmonized_lineage, sep = "||")
  group_factor <- factor(selected_meta$group_id, levels = unique(selected_meta$group_id))
  aggregation <- sparse.model.matrix(~ 0 + group_factor)
  colnames(aggregation) <- levels(group_factor)
  aggregated <- counts[, rownames(selected_meta), drop = FALSE] %*% aggregation
  group_parts <- strsplit(colnames(aggregated), "\\|\\|")
  aggregated_meta <- data.frame(
    pseudobulk_id = colnames(aggregated),
    sample_id = vapply(group_parts, `[[`, character(1), 1L),
    harmonized_lineage = vapply(group_parts, `[[`, character(1), 2L),
    cell_count = as.integer(tapply(as.numeric(cell_weights[keep]), group_factor, sum)[colnames(aggregated)]),
    stringsAsFactors = FALSE
  )
  aggregated_meta <- merge(aggregated_meta, sample_info, by = "sample_id", all.x = TRUE, sort = FALSE)
  aggregated_meta <- aggregated_meta[match(colnames(aggregated), aggregated_meta$pseudobulk_id), , drop = FALSE]
  aggregated_meta$total_cells <- as.integer(total_cells[aggregated_meta$sample_id])
  list(counts = as(aggregated, "dgCMatrix"), meta = aggregated_meta)
}

discovery_count_cache <- file.path(count_root, "discovery_harmonized_raw_counts.rds")
discovery_meta_cache <- file.path(count_root, "discovery_harmonized_metadata.csv")
discovery_t_audit_cache <- file.path(result, "discovery_T_lineage_exclusion_audit.csv")
reuse_discovery <- identical(Sys.getenv("BD_REUSE_STAGE02_DISCOVERY", "0"), "1") &&
  all(file.exists(c(discovery_count_cache, discovery_meta_cache, discovery_t_audit_cache)))
if (reuse_discovery) {
  message("Reusing recovery-only discovery pseudobulk cache")
  discovery_pb <- list(
    counts = as(readRDS(discovery_count_cache), "dgCMatrix"),
    meta = read_csv_path(discovery_meta_cache)
  )
} else {
  message("Loading frozen discovery object")
  discovery_object <- readRDS(discovery_object_path)
  if (!inherits(discovery_object, "Seurat")) stop("Discovery object is not a Seurat object")
  DefaultAssay(discovery_object) <- "RNA"
  discovery_object <- tryCatch(JoinLayers(discovery_object, assay = "RNA"), error = function(e) discovery_object)
  discovery_counts <- tryCatch(LayerData(discovery_object, assay = "RNA", layer = "counts"), error = function(e) NULL)
  if (is.null(discovery_counts)) discovery_counts <- GetAssayData(discovery_object, assay = "RNA", slot = "counts")
  discovery_counts <- as(discovery_counts, "dgCMatrix")
  discovery_meta <- discovery_object[[]]
  required_discovery_meta <- c("sample_id", "severity", "cell_type_broad", "cell_type_subtype")
  missing_discovery_meta <- setdiff(required_discovery_meta, names(discovery_meta))
  if (length(missing_discovery_meta)) stop("Discovery metadata missing: ", paste(missing_discovery_meta, collapse = ", "))
  if (!identical(colnames(discovery_counts), rownames(discovery_meta))) {
    discovery_meta <- discovery_meta[colnames(discovery_counts), , drop = FALSE]
  }
  if (!identical(colnames(discovery_counts), rownames(discovery_meta))) stop("Discovery counts and metadata are misaligned")
  discovery_group <- unname(discovery_broad_map[as.character(discovery_meta$cell_type_broad)])
  t_cells <- discovery_group == "T_lineage" & !is.na(discovery_group)
  discovery_group[t_cells & as.character(discovery_meta$cell_type_subtype) %in% forbidden_t_subtypes] <- NA_character_
  discovery_sample_info <- unique(discovery_meta[, c("sample_id", "severity"), drop = FALSE])
  if (anyDuplicated(discovery_sample_info$sample_id)) stop("Discovery sample metadata is not unique")
  discovery_sample_info$sample_id <- as.character(discovery_sample_info$sample_id)
  discovery_sample_info$severity <- factor(as.character(discovery_sample_info$severity), levels = c("mild", "severe"))
  discovery_total_cells <- table(as.character(discovery_meta$sample_id))
  discovery_pb <- aggregate_cells(
    discovery_counts,
    discovery_meta,
    discovery_group,
    discovery_sample_info,
    discovery_total_cells
  )
  saveRDS(discovery_pb$counts, discovery_count_cache, compress = "gzip")
  write.csv(discovery_pb$meta, discovery_meta_cache, row.names = FALSE)

  discovery_t_audit <- data.frame(
    sample_id = sort(unique(as.character(discovery_meta$sample_id))),
    stringsAsFactors = FALSE
  )
  discovery_t_audit$total_T_NK <- as.integer(table(as.character(discovery_meta$sample_id[discovery_meta$cell_type_broad == "T_NK"]))[discovery_t_audit$sample_id])
  discovery_t_audit$total_T_NK[is.na(discovery_t_audit$total_T_NK)] <- 0L
  kept_t <- discovery_meta$cell_type_broad == "T_NK" & !as.character(discovery_meta$cell_type_subtype) %in% forbidden_t_subtypes
  discovery_t_audit$kept_T_lineage <- as.integer(table(as.character(discovery_meta$sample_id[kept_t]))[discovery_t_audit$sample_id])
  discovery_t_audit$kept_T_lineage[is.na(discovery_t_audit$kept_T_lineage)] <- 0L
  discovery_t_audit$excluded_NK_Cycling_Unresolved <- discovery_t_audit$total_T_NK - discovery_t_audit$kept_T_lineage
  write_csv(discovery_t_audit, "discovery_T_lineage_exclusion_audit.csv")
  rm(discovery_object, discovery_counts, discovery_meta, discovery_group)
  invisible(gc(verbose = FALSE))
}

message("Loading validation pseudobulk raw counts")
validation_counts <- readRDS(validation_count_path)
validation_counts <- as(validation_counts, "dgCMatrix")
validation_meta <- read_csv_path(validation_meta_path)
if (!all(c("pseudobulk_id", "sample_id", "celltype", "cell_count", "severity", "allergy_group", "plastic_group", "batch_wave") %in% names(validation_meta))) {
  stop("Validation pseudobulk metadata schema is incomplete")
}
validation_meta <- validation_meta[match(colnames(validation_counts), validation_meta$pseudobulk_id), , drop = FALSE]
if (anyNA(validation_meta$pseudobulk_id)) stop("Validation counts and metadata are misaligned")
validation_sample_info <- unique(validation_meta[, c("sample_id", "severity", "allergy_group", "plastic_group", "batch_wave"), drop = FALSE])
if (anyDuplicated(validation_sample_info$sample_id)) stop("Validation sample metadata is not unique")
validation_sample_info$severity <- factor(
  ifelse(validation_sample_info$severity == "\u91cd\u75c7", "severe", "mild"),
  levels = c("mild", "severe")
)
validation_sample_info$allergy_group <- factor(
  ifelse(validation_sample_info$allergy_group == "\u662f", "yes", "no"),
  levels = c("no", "yes")
)
validation_sample_info$plastic_group <- factor(
  ifelse(validation_sample_info$plastic_group == "\u662f", "yes", "no"),
  levels = c("no", "yes")
)
validation_sample_info$batch_wave <- factor(validation_sample_info$batch_wave)
validation_total_cells <- aggregate(cell_count ~ sample_id, validation_meta, sum)
validation_total_cells <- setNames(validation_total_cells$cell_count, validation_total_cells$sample_id)

validation_group <- rep(NA_character_, nrow(validation_meta))
for (lineage in names(validation_label_map)) {
  validation_group[validation_meta$celltype %in% validation_label_map[[lineage]]] <- lineage
}
rownames(validation_meta) <- validation_meta$pseudobulk_id
validation_pb <- aggregate_cells(
  validation_counts,
  validation_meta,
  validation_group,
  validation_sample_info,
  validation_total_cells,
  cell_weights = validation_meta$cell_count
)
saveRDS(validation_pb$counts, file.path(count_root, "validation_harmonized_raw_counts.rds"), compress = "gzip")
write.csv(validation_pb$meta, file.path(count_root, "validation_harmonized_metadata.csv"), row.names = FALSE)

mac_source <- validation_meta[validation_meta$celltype %in% c("Macrophage", "Monocyte"), c("sample_id", "celltype", "cell_count"), drop = FALSE]
mac_source_wide <- reshape(mac_source, idvar = "sample_id", timevar = "celltype", direction = "wide")
names(mac_source_wide) <- sub("^cell_count\\.", "", names(mac_source_wide))
mac_aggregated <- validation_pb$meta[validation_pb$meta$harmonized_lineage == "Macrophage_Monocyte", c("sample_id", "cell_count"), drop = FALSE]
names(mac_aggregated)[[2]] <- "aggregated_cell_count"
mac_audit <- merge(mac_source_wide, mac_aggregated, by = "sample_id", all = TRUE)
for (column in c("Macrophage", "Monocyte")) if (!column %in% names(mac_audit)) mac_audit[[column]] <- 0L
mac_audit$source_sum <- rowSums(mac_audit[, c("Macrophage", "Monocyte"), drop = FALSE], na.rm = TRUE)
mac_audit$exact_sum <- mac_audit$source_sum == mac_audit$aggregated_cell_count
write_csv(mac_audit, "validation_macrophage_monocyte_sum_audit.csv")

t_markers <- intersect(c("CD3D", "CD3E", "CD3G", "TRAC", "LCK"), rownames(validation_counts))
if (length(t_markers) < 4L) stop("Too few T-lineage marker genes in validation counts")
validation_logcpm <- cpm(DGEList(validation_counts), log = TRUE, prior.count = 1)
t_score <- colMeans(validation_logcpm[t_markers, , drop = FALSE])
t_detail <- data.frame(
  pseudobulk_id = colnames(validation_counts),
  sample_id = validation_meta$sample_id,
  celltype = validation_meta$celltype,
  t_marker_score = as.numeric(t_score),
  cd3_trac_raw_count = as.numeric(colSums(validation_counts[intersect(c("CD3D", "CD3E", "TRAC"), rownames(validation_counts)), , drop = FALSE])),
  stringsAsFactors = FALSE
)
t_columns <- t_detail$celltype == "T_cells"
t_support <- data.frame(
  marker_genes = paste(t_markers, collapse = ";"),
  t_pseudobulks = sum(t_columns),
  cd3_trac_positive_fraction = mean(t_detail$cd3_trac_raw_count[t_columns] > 0),
  median_t_marker_score_T_cells = median(t_detail$t_marker_score[t_columns]),
  median_t_marker_score_non_T = median(t_detail$t_marker_score[!t_columns]),
  median_score_delta = median(t_detail$t_marker_score[t_columns]) - median(t_detail$t_marker_score[!t_columns]),
  support_status = if (mean(t_detail$cd3_trac_raw_count[t_columns] > 0) >= 0.9 && median(t_detail$t_marker_score[t_columns]) - median(t_detail$t_marker_score[!t_columns]) > 1) "PASS" else "FAIL",
  stringsAsFactors = FALSE
)
write_csv(t_detail, "validation_T_lineage_marker_detail.csv")
write_csv(t_support, "validation_T_lineage_marker_support.csv")

prepare_lineage <- function(pb, lineage, adjusted = FALSE) {
  columns <- which(pb$meta$harmonized_lineage == lineage & pb$meta$cell_count >= 50L)
  meta <- pb$meta[columns, , drop = FALSE]
  counts <- pb$counts[, columns, drop = FALSE]
  meta$sample_id <- as.character(meta$sample_id)
  meta$severity <- factor(as.character(meta$severity), levels = c("mild", "severe"))
  if (adjusted) {
    meta$batch_wave <- factor(meta$batch_wave)
    meta$allergy_group <- factor(meta$allergy_group, levels = c("no", "yes"))
    meta$plastic_group <- factor(meta$plastic_group, levels = c("no", "yes"))
  }
  list(counts = counts, meta = meta)
}

run_de_model <- function(counts, meta, cohort, lineage, model_id, adjusted = FALSE) {
  group_counts <- table(factor(meta$severity, levels = c("mild", "severe")))
  if (any(group_counts < 3L)) stop(cohort, " ", lineage, " has fewer than three eligible patients per severity group")
  design <- if (adjusted) {
    model.matrix(~ batch_wave + severity + allergy_group + plastic_group, data = meta)
  } else {
    model.matrix(~ severity, data = meta)
  }
  if (qr(design)$rank != ncol(design)) stop(cohort, " ", lineage, " model is not full rank: ", model_id)
  coefficient <- "severitysevere"
  if (!coefficient %in% colnames(design)) stop("Severity coefficient is absent: ", cohort, " ", lineage)

  y <- DGEList(counts = counts, samples = meta, group = meta$severity)
  keep <- filterByExpr(y, group = meta$severity, min.count = 10, min.total.count = 15)
  if (sum(keep) < 100L) stop(cohort, " ", lineage, " has too few genes after filtering")
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")
  y <- estimateDisp(y, design, robust = TRUE)
  ql_fit <- glmQLFit(y, design, robust = TRUE)
  ql_test <- glmQLFTest(ql_fit, coef = coefficient)
  edge <- topTags(ql_test, n = Inf, sort.by = "none")$table
  edge$gene <- rownames(edge)

  voomed <- voom(y, design, plot = FALSE)
  voom_fit <- eBayes(lmFit(voomed, design), robust = TRUE)
  coefficient_index <- match(coefficient, colnames(voom_fit$coefficients))
  effect <- voom_fit$coefficients[, coefficient_index]
  effect_se <- voom_fit$stdev.unscaled[, coefficient_index] * sqrt(voom_fit$s2.post)
  voom_frame <- data.frame(
    gene = rownames(voom_fit$coefficients),
    effect_log2fc = as.numeric(effect),
    effect_se = as.numeric(effect_se),
    stringsAsFactors = FALSE
  )
  voom_frame$effect_ci95_low <- voom_frame$effect_log2fc - 1.96 * voom_frame$effect_se
  voom_frame$effect_ci95_high <- voom_frame$effect_log2fc + 1.96 * voom_frame$effect_se
  voom_frame <- voom_frame[match(edge$gene, voom_frame$gene), , drop = FALSE]

  result_frame <- data.frame(
    cohort = cohort,
    harmonized_lineage = lineage,
    model_id = model_id,
    contrast = "severe_vs_mild",
    gene = edge$gene,
    edgeR_logFC = edge$logFC,
    effect_log2fc = voom_frame$effect_log2fc,
    effect_se = voom_frame$effect_se,
    effect_ci95_low = voom_frame$effect_ci95_low,
    effect_ci95_high = voom_frame$effect_ci95_high,
    logCPM = edge$logCPM,
    F = edge$F,
    PValue = edge$PValue,
    FDR = edge$FDR,
    n_mild_patients = as.integer(group_counts[["mild"]]),
    n_severe_patients = as.integer(group_counts[["severe"]]),
    min_cells_per_pseudobulk = min(meta$cell_count),
    design_rank = qr(design)$rank,
    design_columns = ncol(design),
    effect_method = "limma_voom_empirical_bayes",
    test_method = "edgeR_robust_QLF",
    stringsAsFactors = FALSE
  )
  result_frame <- result_frame[order(result_frame$FDR, -abs(result_frame$effect_log2fc)), , drop = FALSE]
  summary_frame <- data.frame(
    cohort = cohort,
    harmonized_lineage = lineage,
    model_id = model_id,
    samples_used = nrow(meta),
    n_mild_patients = as.integer(group_counts[["mild"]]),
    n_severe_patients = as.integer(group_counts[["severe"]]),
    min_cells_per_pseudobulk = min(meta$cell_count),
    genes_tested = nrow(result_frame),
    fdr_lt_0_05 = sum(result_frame$FDR < 0.05),
    fdr_lt_0_05_abs_effect_ge_0_5 = sum(result_frame$FDR < 0.05 & abs(result_frame$effect_log2fc) >= 0.5),
    finite_positive_effect_se = sum(is.finite(result_frame$effect_se) & result_frame$effect_se > 0),
    design_rank = qr(design)$rank,
    design_columns = ncol(design),
    status = "COMPLETED",
    stringsAsFactors = FALSE
  )
  list(result = result_frame, summary = summary_frame)
}

primary_results <- list()
primary_summaries <- list()
adjusted_results <- list()
adjusted_summaries <- list()
index <- 0L
adjusted_index <- 0L
for (lineage in lineages) {
  index <- index + 1L
  discovery_current <- prepare_lineage(discovery_pb, lineage)
  primary_results[[index]] <- run_de_model(
    discovery_current$counts, discovery_current$meta, "discovery_10", lineage, "harmonized_severity_only"
  )
  index <- index + 1L
  validation_current <- prepare_lineage(validation_pb, lineage)
  primary_results[[index]] <- run_de_model(
    validation_current$counts, validation_current$meta, "validation_23", lineage, "harmonized_severity_only"
  )
  adjusted_index <- adjusted_index + 1L
  validation_adjusted <- prepare_lineage(validation_pb, lineage, adjusted = TRUE)
  adjusted_results[[adjusted_index]] <- run_de_model(
    validation_adjusted$counts,
    validation_adjusted$meta,
    "validation_23",
    lineage,
    "validation_covariate_adjusted",
    adjusted = TRUE
  )
}
primary_de <- do.call(rbind, lapply(primary_results, `[[`, "result"))
primary_summary <- do.call(rbind, lapply(primary_results, `[[`, "summary"))
adjusted_de <- do.call(rbind, lapply(adjusted_results, `[[`, "result"))
adjusted_summary <- do.call(rbind, lapply(adjusted_results, `[[`, "summary"))
write.csv(primary_de, gzfile(file.path(result, "harmonized_primary_de.csv.gz")), row.names = FALSE)
write_csv(primary_summary, "harmonized_primary_de_summary.csv")
write.csv(adjusted_de, gzfile(file.path(result, "validation_adjusted_de.csv.gz")), row.names = FALSE)
write_csv(adjusted_summary, "validation_adjusted_de_summary.csv")

composition_columns <- c(
  "sample_id", "severity", "harmonized_lineage", "cell_count", "total_cells"
)
composition <- rbind(
  transform(discovery_pb$meta[, composition_columns, drop = FALSE], cohort = "discovery_10"),
  transform(validation_pb$meta[, composition_columns, drop = FALSE], cohort = "validation_23")
)
composition$proportion <- composition$cell_count / composition$total_cells
composition <- composition[, c("cohort", "sample_id", "severity", "harmonized_lineage", "cell_count", "total_cells", "proportion")]
write_csv(composition, "harmonized_composition_by_sample.csv")

bootstrap_difference <- function(mild, severe, iterations = 5000L) {
  differences <- replicate(iterations, mean(sample(severe, replace = TRUE)) - mean(sample(mild, replace = TRUE)))
  unname(quantile(differences, c(0.025, 0.975), na.rm = TRUE))
}
composition_results <- do.call(rbind, lapply(split(composition, list(composition$cohort, composition$harmonized_lineage), drop = TRUE), function(frame) {
  mild <- frame$proportion[frame$severity == "mild"]
  severe <- frame$proportion[frame$severity == "severe"]
  interval <- bootstrap_difference(mild, severe)
  test <- suppressWarnings(wilcox.test(severe, mild, exact = FALSE))
  effect_se <- sqrt(var(severe) / length(severe) + var(mild) / length(mild))
  data.frame(
    cohort = frame$cohort[[1]],
    harmonized_lineage = frame$harmonized_lineage[[1]],
    contrast = "severe_minus_mild_proportion",
    n_mild_patients = length(mild),
    n_severe_patients = length(severe),
    mild_mean = mean(mild),
    severe_mean = mean(severe),
    mean_difference = mean(severe) - mean(mild),
    effect_se = effect_se,
    ci95_low = mean(severe) - mean(mild) - 1.96 * effect_se,
    ci95_high = mean(severe) - mean(mild) + 1.96 * effect_se,
    bootstrap_ci_low = interval[[1]],
    bootstrap_ci_high = interval[[2]],
    wilcoxon_p_value = test$p.value,
    stringsAsFactors = FALSE
  )
}))
composition_results$fdr_within_cohort <- ave(
  composition_results$wilcoxon_p_value,
  composition_results$cohort,
  FUN = function(x) p.adjust(x, method = "BH")
)
write_csv(composition_results, "harmonized_composition_comparison.csv")

input_paths <- c(
  discovery_object_path,
  file.path(freeze_root, "discovery_v1_freeze_record.csv"),
  validation_count_path,
  validation_meta_path
)
input_info <- file.info(input_paths)
input_inventory <- data.frame(
  path = input_paths,
  exists = file.exists(input_paths),
  bytes = as.numeric(input_info$size),
  modified = format(input_info$mtime, "%Y-%m-%dT%H:%M:%S%z"),
  integrity_reference = c(
    "discovery_v1_manifest",
    "stage_01_freeze_validation",
    "stage_00_input_inventory",
    "stage_00_input_inventory"
  ),
  stringsAsFactors = FALSE
)
write_csv(input_inventory, "stage_02_input_inventory.csv")

report <- c(
  "# Harmonized discovery-validation analysis",
  "",
  paste0("- Generated: ", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
  "- Statistical replicate: patient sample.",
  "- Cohorts are independent and are never patient-matched.",
  "- Harmonized lineages: B, macrophage/monocyte, neutrophil, and T.",
  "- Discovery T lineage excludes NK, Cycling_T_NK, and Unresolved cells before aggregation.",
  "- Validation macrophage and monocyte raw counts are summed by patient before modeling.",
  "- Primary comparable design in both cohorts: ~ severity.",
  "- Validation sensitivity design: ~ batch_wave + severity + allergy_group + plastic_group.",
  "- Differential-expression test: edgeR robust quasi-likelihood.",
  "- Meta-compatible effect and standard error: limma-voom empirical Bayes coefficient and moderated standard error.",
  "- Pseudobulks require at least 50 cells; each severity group requires at least three patients.",
  "- Gene filtering: edgeR filterByExpr with min.count=10 and min.total.count=15."
)
writeLines(report, file.path(result, "harmonized_validation_report.md"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(result, "sessionInfo.txt"), useBytes = TRUE)
message("Harmonized validation analysis completed")
