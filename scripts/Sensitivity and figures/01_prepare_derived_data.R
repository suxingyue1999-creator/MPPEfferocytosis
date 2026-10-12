#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
})

set.seed(20260728)

workspace <- Sys.getenv("MPP_WORKSPACE", "inputs/workspace")
package_dir <- file.path(workspace, "工作文件", "论文主图_2026-07-28")
data_dir <- file.path(package_dir, "data")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

discovery_root <- file.path(
  workspace,
  "BD_Rhapsody_setup/mac_portable_full/bd_rhapsody_10/downstream",
  "10samples_scrna_mex_v4/results"
)
validation_root <- file.path(
  workspace,
  "工作文件/23例单细胞重新分析/results"
)
cross_root <- file.path(
  workspace,
  "BD_Rhapsody_setup/figure_work_package/cross_cohort_results"
)

required_files <- c(
  discovery_metadata = file.path(discovery_root, "00_audit/10samples_metadata_locked.csv"),
  discovery_qc = file.path(discovery_root, "01_qc/10samples_qc_summary.csv"),
  discovery_umap = file.path(discovery_root, "02_integrate/cell_embeddings.csv.gz"),
  discovery_annotations = file.path(discovery_root, "03_annotate/cell_annotations.csv.gz"),
  discovery_myeloid_umap = file.path(discovery_root, "03_annotate/lineage_recluster_embeddings.csv.gz"),
  validation_qc = file.path(validation_root, "01_qc/qc_summary_by_sample.csv"),
  validation_atlas_rds = file.path(
    validation_root,
    "02_rna_reanalysis/23samples_RNA_reclustered_annotated.rds"
  ),
  validation_myeloid_rds = file.path(
    validation_root,
    "02_rna_reanalysis/subclusters/Macrophage_Monocyte",
    "Macrophage_Monocyte_subclustered.rds"
  ),
  discovery_pb_counts = file.path(
    cross_root,
    "02_harmonized_validation/pseudobulk_counts/discovery_harmonized_raw_counts.rds"
  ),
  discovery_pb_metadata = file.path(
    cross_root,
    "02_harmonized_validation/pseudobulk_counts/discovery_harmonized_metadata.csv"
  ),
  validation_pb_counts = file.path(
    cross_root,
    "02_harmonized_validation/pseudobulk_counts/validation_harmonized_raw_counts.rds"
  ),
  validation_pb_metadata = file.path(
    cross_root,
    "02_harmonized_validation/pseudobulk_counts/validation_harmonized_metadata.csv"
  ),
  pathway_scores = file.path(
    cross_root,
    "04_replicated_mechanism/replicated_pathway_scores_by_sample.csv"
  ),
  pathway_models = file.path(
    cross_root,
    "04_replicated_mechanism/replicated_pathway_score_models.csv"
  ),
  pathway_meta = file.path(
    cross_root,
    "04_replicated_mechanism/replicated_pathway_score_meta_analysis.csv"
  ),
  myeloid_state_composition = file.path(
    cross_root,
    "04_replicated_mechanism/harmonized_macrophage_state_composition.csv"
  ),
  myeloid_state_models = file.path(
    cross_root,
    "04_replicated_mechanism/harmonized_macrophage_state_models.csv"
  ),
  myeloid_state_meta = file.path(
    cross_root,
    "04_replicated_mechanism/harmonized_macrophage_state_meta_analysis.csv"
  ),
  discovery_qc_sensitivity = file.path(
    discovery_root,
    "04_compare/pseudobulk_sensitivity_concordance.csv"
  ),
  validation_covariate_concordance = file.path(
    cross_root,
    "03_meta_pathway/validation_adjustment_concordance.csv"
  ),
  composition_meta = file.path(
    cross_root,
    "03_meta_pathway/two_cohort_composition_meta_analysis.csv"
  )
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Missing required input files:\n", paste(missing_files, collapse = "\n"))
}

message("Preparing discovery atlas coordinates and conservative lineage labels")
discovery_umap <- fread(required_files[["discovery_umap"]])
discovery_annotations <- fread(required_files[["discovery_annotations"]])
if (!identical(discovery_umap$cell_id, discovery_annotations$cell_id)) {
  stop("Discovery UMAP and annotation cell order differ")
}

discovery_atlas <- discovery_umap[, .(cell_id, umap_1, umap_2, sample_id, severity)]
discovery_atlas[, cluster := as.character(discovery_annotations$seurat_clusters)]
discovery_atlas[, source_label := discovery_annotations$cell_type_broad]
discovery_atlas[, lineage := fcase(
  source_label == "B_Plasma", "B lineage",
  source_label == "Macrophage_Monocyte", "Macrophage/monocyte",
  source_label == "Neutrophil", "Neutrophil",
  source_label == "T_NK", "T lineage",
  default = "Other"
)]

# Marker review identified four rare technical clusters that were mislabeled by
# the automated broad annotation. Restrict the atlas to conservative lineages.
discovery_atlas[cluster == "14", lineage := "Macrophage/monocyte"]
discovery_atlas[cluster %in% c("17", "23"), lineage := "T lineage"]
discovery_atlas[cluster == "25", lineage := "B lineage"]
discovery_atlas[, cohort := "Discovery"]
fwrite(
  discovery_atlas[, .(cell_id, umap_1, umap_2, sample_id, severity, lineage, cohort)],
  file.path(data_dir, "discovery_atlas_umap.csv.gz")
)

message("Extracting validation atlas coordinates from the current RNA object")
validation_object <- readRDS(required_files[["validation_atlas_rds"]])
validation_reduction <- "umap.rpca"
if (!validation_reduction %in% Reductions(validation_object)) {
  stop("Expected reduction not found: ", validation_reduction)
}
validation_coordinates <- as.data.table(
  Embeddings(validation_object, validation_reduction),
  keep.rownames = "cell_id"
)
setnames(validation_coordinates, 2:3, c("umap_1", "umap_2"))
validation_metadata <- as.data.table(
  validation_object@meta.data[, c("sample_id", "severity", "celltype_broad_rechecked")],
  keep.rownames = "cell_id"
)
validation_atlas <- merge(
  validation_coordinates,
  validation_metadata,
  by = "cell_id",
  sort = FALSE
)
if (nrow(validation_atlas) != ncol(validation_object)) {
  stop("Validation atlas coordinate/metadata join lost cells")
}
validation_atlas[, severity := fcase(
  severity %chin% c("轻症", "mild"), "mild",
  severity %chin% c("重症", "severe"), "severe",
  default = NA_character_
)]
validation_atlas[, lineage := fcase(
  celltype_broad_rechecked == "B_cell", "B lineage",
  celltype_broad_rechecked %chin% c("Macrophage", "Monocyte"), "Macrophage/monocyte",
  celltype_broad_rechecked == "Neutrophils", "Neutrophil",
  celltype_broad_rechecked == "T_cells", "T lineage",
  default = "Other"
)]
validation_atlas[, cohort := "Validation"]
fwrite(
  validation_atlas[, .(cell_id, umap_1, umap_2, sample_id, severity, lineage, cohort)],
  file.path(data_dir, "validation_atlas_umap.csv.gz")
)
rm(validation_object, validation_coordinates, validation_metadata, validation_atlas)
invisible(gc())

message("Preparing discovery myeloid-state coordinates")
discovery_myeloid <- fread(required_files[["discovery_myeloid_umap"]])[
  lineage == "Macrophage_Monocyte"
]
discovery_myeloid <- merge(
  discovery_myeloid,
  discovery_annotations[, .(cell_id, sample_id, severity)],
  by = "cell_id",
  sort = FALSE
)
discovery_myeloid[, state := fcase(
  cell_type_subtype == "Alveolar_Macrophage", "Alveolar macrophages",
  cell_type_subtype == "Cycling_Myeloid", "Cycling myeloid",
  default = "Other myeloid"
)]
discovery_myeloid[, cohort := "Discovery"]
fwrite(
  discovery_myeloid[, .(cell_id, umap_1, umap_2, sample_id, severity, state, cohort)],
  file.path(data_dir, "discovery_myeloid_umap.csv.gz")
)

message("Extracting validation myeloid-state coordinates")
validation_myeloid_object <- readRDS(required_files[["validation_myeloid_rds"]])
validation_myeloid_reduction <- "umap.lineage"
if (!validation_myeloid_reduction %in% Reductions(validation_myeloid_object)) {
  stop("Expected reduction not found: ", validation_myeloid_reduction)
}
validation_myeloid_coordinates <- as.data.table(
  Embeddings(validation_myeloid_object, validation_myeloid_reduction),
  keep.rownames = "cell_id"
)
setnames(validation_myeloid_coordinates, 2:3, c("umap_1", "umap_2"))
validation_myeloid_metadata <- as.data.table(
  validation_myeloid_object@meta.data[, c("sample_id", "severity", "subtype_label")],
  keep.rownames = "cell_id"
)
validation_myeloid <- merge(
  validation_myeloid_coordinates,
  validation_myeloid_metadata,
  by = "cell_id",
  sort = FALSE
)
if (nrow(validation_myeloid) != ncol(validation_myeloid_object)) {
  stop("Validation myeloid coordinate/metadata join lost cells")
}
validation_myeloid[, severity := fcase(
  severity %chin% c("轻症", "mild"), "mild",
  severity %chin% c("重症", "severe"), "severe",
  default = NA_character_
)]
validation_myeloid[, state := fcase(
  subtype_label == "FABP4_alveolar", "Alveolar macrophages",
  subtype_label == "proliferating", "Cycling myeloid",
  default = "Other myeloid"
)]
validation_myeloid[, cohort := "Validation"]
fwrite(
  validation_myeloid[, .(cell_id, umap_1, umap_2, sample_id, severity, state, cohort)],
  file.path(data_dir, "validation_myeloid_umap.csv.gz")
)
rm(
  validation_myeloid_object,
  validation_myeloid_coordinates,
  validation_myeloid_metadata,
  validation_myeloid
)
invisible(gc())

message("Combining patient-level QC retention")
discovery_qc <- fread(required_files[["discovery_qc"]])[
  , .(
    cohort = "Discovery",
    sample_id,
    severity,
    cells_input,
    cells_retained = cells_qc_keep,
    retention_rate
  )
]
validation_qc <- fread(required_files[["validation_qc"]])[
  , .(
    cohort = "Validation",
    sample_id,
    severity = fcase(
      severity %chin% c("轻症", "mild"), "mild",
      severity %chin% c("重症", "severe"), "severe",
      default = NA_character_
    ),
    cells_input = cells_raw_matrix,
    cells_retained = cells_qc_pass,
    retention_rate = qc_pass_rate
  )
]
qc_retention <- rbindlist(list(discovery_qc, validation_qc), use.names = TRUE)
fwrite(qc_retention, file.path(data_dir, "qc_retention.csv"))

message("Computing canonical marker expression from patient pseudobulks")
marker_genes <- c(
  "MS4A1", "CD79A", "CD74",
  "LST1", "FCER1G", "CTSS",
  "FCGR3B", "CSF3R", "S100A8",
  "CD3D", "CD3E", "TRBC1"
)

marker_groups <- data.table(
  gene = marker_genes,
  marker_lineage = rep(
    c("B lineage", "Macrophage/monocyte", "Neutrophil", "T lineage"),
    each = 3
  )
)

summarize_marker_expression <- function(counts_path, metadata_path, cohort_name) {
  counts <- readRDS(counts_path)
  metadata <- fread(metadata_path)
  if (!all(marker_genes %in% rownames(counts))) {
    stop("Missing canonical genes in ", cohort_name, ": ",
         paste(setdiff(marker_genes, rownames(counts)), collapse = ", "))
  }
  metadata <- metadata[match(colnames(counts), pseudobulk_id)]
  if (anyNA(metadata$pseudobulk_id)) {
    stop("Pseudobulk metadata does not match count columns for ", cohort_name)
  }
  library_sizes <- Matrix::colSums(counts)
  marker_counts <- as.matrix(counts[marker_genes, , drop = FALSE])
  marker_log_cpm <- log2(sweep(marker_counts, 2, library_sizes, "/") * 1e6 + 1)
  marker_long <- as.data.table(as.table(marker_log_cpm))
  setnames(marker_long, c("gene", "pseudobulk_id", "log2_cpm"))
  marker_long[, gene := as.character(gene)]
  marker_long[, pseudobulk_id := as.character(pseudobulk_id)]
  marker_long <- merge(
    marker_long,
    metadata[, .(pseudobulk_id, harmonized_lineage)],
    by = "pseudobulk_id",
    sort = FALSE
  )
  result <- marker_long[
    , .(
      mean_log2_cpm = mean(log2_cpm),
      sd_log2_cpm = sd(log2_cpm),
      patients = .N
    ),
    by = .(gene, harmonized_lineage)
  ]
  result[, cohort := cohort_name]
  result[]
}

marker_expression <- rbindlist(list(
  summarize_marker_expression(
    required_files[["discovery_pb_counts"]],
    required_files[["discovery_pb_metadata"]],
    "Discovery"
  ),
  summarize_marker_expression(
    required_files[["validation_pb_counts"]],
    required_files[["validation_pb_metadata"]],
    "Validation"
  )
))
marker_expression <- merge(marker_expression, marker_groups, by = "gene")
lineage_labels <- c(
  B_lineage = "B lineage",
  Macrophage_Monocyte = "Macrophage/monocyte",
  Neutrophil = "Neutrophil",
  T_lineage = "T lineage"
)
marker_expression[, lineage := unname(lineage_labels[harmonized_lineage])]
marker_expression[
  , expression_z := as.numeric(scale(mean_log2_cpm)),
  by = gene
]
fwrite(marker_expression, file.path(data_dir, "canonical_marker_pseudobulk_expression.csv"))

message("Recomputing pathway leave-one-patient-out sensitivity")
pathway_scores <- fread(required_files[["pathway_scores"]])
pathway_models <- fread(required_files[["pathway_models"]])
pathway_scores[, severity := factor(severity, levels = c("mild", "severe"))]
pathway_scores[, batch_wave := factor(batch_wave)]
pathway_scores[, allergy_group := factor(allergy_group, levels = c("no", "yes"))]
pathway_scores[, plastic_group := factor(plastic_group, levels = c("no", "yes"))]

fit_pathway_model <- function(data, model_key) {
  model_formula <- if (model_key == "validation_adjusted") {
    standardized_score ~ batch_wave + severity + allergy_group + plastic_group
  } else {
    standardized_score ~ severity
  }
  fit <- lm(model_formula, data = data)
  coefficient <- unname(coef(fit)["severitysevere"])
  if (!is.finite(coefficient)) {
    stop("Non-estimable severity coefficient in ", model_key)
  }
  coefficient
}

loo_specs <- data.table(
  cohort = c("discovery_10", "validation_23", "validation_23"),
  model_key = c("discovery_primary", "validation_primary", "validation_adjusted"),
  model_id = c("severity_only", "severity_only", "covariate_adjusted")
)

loo_records <- vector("list", nrow(loo_specs) * uniqueN(pathway_scores$go_id))
record_index <- 1L
for (spec_index in seq_len(nrow(loo_specs))) {
  spec <- loo_specs[spec_index]
  cohort_scores <- pathway_scores[cohort == spec$cohort]
  for (pathway_id in unique(cohort_scores$go_id)) {
    pathway_data <- cohort_scores[go_id == pathway_id]
    full_effect <- fit_pathway_model(pathway_data, spec$model_key)
    published_effect <- pathway_models[
      cohort == spec$cohort &
        model_id == spec$model_id &
        go_id == pathway_id &
        coefficient == "severitysevere",
      estimate
    ]
    if (length(published_effect) != 1L || abs(full_effect - published_effect) > 1e-8) {
      stop("LOO full-fit coefficient does not match published model for ", spec$model_key,
           " / ", pathway_id)
    }
    omissions <- lapply(unique(pathway_data$sample_id), function(omitted_id) {
      effect <- fit_pathway_model(
        pathway_data[sample_id != omitted_id],
        spec$model_key
      )
      data.table(
        cohort = spec$cohort,
        model_key = spec$model_key,
        go_id = pathway_id,
        term = pathway_data$term[[1]],
        omitted_sample = omitted_id,
        full_effect = full_effect,
        loo_effect = effect
      )
    })
    loo_records[[record_index]] <- rbindlist(omissions)
    record_index <- record_index + 1L
  }
}
pathway_loo <- rbindlist(loo_records, use.names = TRUE)
pathway_loo_summary <- pathway_loo[
  , .(
    full_effect = unique(full_effect),
    omissions = .N,
    min_loo_effect = min(loo_effect),
    max_loo_effect = max(loo_effect),
    sign_concordance = mean(sign(loo_effect) == sign(full_effect)),
    all_negative = all(loo_effect < 0)
  ),
  by = .(cohort, model_key, go_id, term)
]
fwrite(pathway_loo, file.path(data_dir, "pathway_leave_one_patient_out.csv.gz"))
fwrite(pathway_loo_summary, file.path(data_dir, "pathway_loo_summary.csv"))

message("Writing source manifest")
source_manifest <- data.table(
  source_id = names(required_files),
  path = unname(required_files),
  bytes = unname(file.info(required_files)$size),
  modified = format(file.info(required_files)$mtime, "%Y-%m-%d %H:%M:%S %z"),
  md5 = unname(tools::md5sum(required_files))
)
fwrite(source_manifest, file.path(package_dir, "source_manifest.csv"))

validation <- data.table(
  check = c(
    "discovery_atlas_cells",
    "validation_atlas_cells",
    "discovery_myeloid_cells",
    "validation_myeloid_cells",
    "qc_patients",
    "marker_rows",
    "marker_missing_values",
    "pathway_loo_models",
    "pathway_loo_rows",
    "pathway_loo_all_negative"
  ),
  observed = c(
    nrow(discovery_atlas),
    nrow(fread(file.path(data_dir, "validation_atlas_umap.csv.gz"))),
    nrow(discovery_myeloid),
    nrow(fread(file.path(data_dir, "validation_myeloid_umap.csv.gz"))),
    nrow(qc_retention),
    nrow(marker_expression),
    sum(is.na(marker_expression$expression_z)),
    nrow(pathway_loo_summary),
    nrow(pathway_loo),
    sum(pathway_loo_summary$all_negative)
  ),
  expected = c(115594, 209299, 23228, 118518, 33, 96, 0, 21, 385, 21)
)
validation[, status := fifelse(observed == expected, "PASS", "FAIL")]
fwrite(validation, file.path(data_dir, "derived_data_validation.csv"))

if (any(validation$status == "FAIL")) {
  stop("Derived-data validation failed; inspect derived_data_validation.csv")
}

message("Derived data preparation completed successfully")
