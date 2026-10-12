#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
})

set.seed(20260721)
options(future.globals.maxSize = 40 * 1024^3)
future::plan("sequential")

workspace <- Sys.getenv("MPP_WORKSPACE", "inputs/workspace")
work_root <- file.path(workspace, "工作文件", "23例单细胞重新分析")
input_path <- file.path(
  work_root,
  "results",
  "01_qc",
  "23samples_current_RNA_qc_counts.rds"
)
out_dir <- file.path(work_root, "results", "02_rna_reanalysis")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_path)) stop("QC object is missing: ", input_path)

message("Loading RNA QC object")
object <- readRDS(input_path)
DefaultAssay(object) <- "RNA"
if (ncol(object) != 209299L) {
  warning("Expected 209,299 QC cells, found ", ncol(object))
}

message("Log-normalizing each sample layer")
object <- NormalizeData(
  object,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = TRUE
)
message("Selecting consensus variable genes")
object <- FindVariableFeatures(
  object,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = TRUE
)
variable_features <- VariableFeatures(object)
writeLines(variable_features, file.path(out_dir, "variable_features_2000.txt"))

message("Scaling variable genes and running PCA")
object <- ScaleData(
  object,
  features = variable_features,
  verbose = TRUE
)
object <- RunPCA(
  object,
  features = variable_features,
  npcs = 50,
  verbose = TRUE
)

pca_stdev <- Stdev(object, reduction = "pca")
pca_variance <- pca_stdev^2 / sum(pca_stdev^2)
pca_table <- data.frame(
  PC = seq_along(pca_stdev),
  stdev = pca_stdev,
  variance_fraction = pca_variance,
  cumulative_variance = cumsum(pca_variance)
)
write.csv(
  pca_table,
  file.path(out_dir, "pca_variance.csv"),
  row.names = FALSE
)
gc(verbose = FALSE)

message("Integrating sample layers with RPCA for clustering/visualization only")
data_layers <- Layers(object[["RNA"]], search = "^data\\.")
reference_samples <- c("V01", "V18", "V07", "V14")
reference_indices <- which(data_layers %in% paste0("data.", reference_samples))
if (length(reference_indices) != length(reference_samples)) {
  stop("Balanced RPCA reference layers were not all found")
}
object <- IntegrateLayers(
  object = object,
  method = RPCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.rpca",
  normalization.method = "LogNormalize",
  dims = 1:30,
  reference = reference_indices,
  k.filter = 100,
  k.weight = 50,
  verbose = TRUE
)

message("Building integrated graph, clusters, and UMAP")
object <- FindNeighbors(
  object,
  reduction = "integrated.rpca",
  dims = 1:30,
  k.param = 30,
  verbose = TRUE
)
object <- FindClusters(
  object,
  resolution = c(0.2, 0.4, 0.6, 0.8),
  algorithm = 4,
  leiden_method = "igraph",
  random.seed = 20260721,
  verbose = TRUE
)
object <- RunUMAP(
  object,
  reduction = "integrated.rpca",
  dims = 1:30,
  n.neighbors = 30,
  min.dist = 0.3,
  metric = "cosine",
  seed.use = 20260721,
  reduction.name = "umap.rpca",
  verbose = TRUE
)
object$cluster_rna <- object$RNA_snn_res.0.6
Idents(object) <- "cluster_rna"

# Integration is finished. Join expression layers so downstream marker and
# pseudobulk code sees one RNA count matrix; integrated values remain only in
# the dimensional reduction.
object <- JoinLayers(object, assay = "RNA")

marker_panels <- list(
  Macrophage = c("C1QA", "C1QB", "C1QC", "APOC1", "FABP4", "MARCO", "PPARG", "MRC1"),
  Monocyte = c("LYZ", "FCN1", "CTSS", "LILRB1", "CTSD", "S100A10", "VCAN"),
  T_cells = c("CD3D", "CD3E", "TRAC", "LCK", "IL7R", "CCL5", "CD247"),
  B_cell = c("CD79A", "CD79B", "MS4A1", "CD37", "CD74", "CD22", "MZB1", "JCHAIN"),
  Neutrophils = c("FCGR3B", "CSF3R", "CXCR2", "FPR1", "S100A8", "S100A9", "MNDA"),
  Epithelial_cells = c("EPCAM", "KRT8", "KRT18", "KRT19", "SCGB1A1", "SFTPC", "KRT5")
)
marker_panels <- lapply(marker_panels, function(genes) intersect(genes, rownames(object)))
if (any(lengths(marker_panels) < 3L)) {
  stop("At least one broad marker panel has fewer than three genes")
}

message("Scoring broad cell-type marker panels")
object <- AddModuleScore(
  object,
  features = marker_panels,
  name = "broad_marker_score_",
  assay = "RNA",
  seed = 20260721,
  search = FALSE
)
score_columns <- paste0("broad_marker_score_", seq_along(marker_panels))
score_matrix <- as.matrix(object[[]][, score_columns, drop = FALSE])
colnames(score_matrix) <- names(marker_panels)
object$marker_broad_call <- colnames(score_matrix)[max.col(score_matrix, ties.method = "first")]

cluster_ids <- levels(factor(object$cluster_rna))
annotation_rows <- lapply(cluster_ids, function(cluster_id) {
  index <- which(object$cluster_rna == cluster_id)
  legacy_labels <- object$previous_celltype_main[index]
  legacy_labels <- legacy_labels[!is.na(legacy_labels) & nzchar(legacy_labels)]
  marker_labels <- object$marker_broad_call[index]
  legacy_table <- sort(table(legacy_labels), decreasing = TRUE)
  marker_table <- sort(table(marker_labels), decreasing = TRUE)
  legacy_majority <- if (length(legacy_table)) names(legacy_table)[1] else NA_character_
  marker_majority <- names(marker_table)[1]
  legacy_purity <- if (length(legacy_table)) as.numeric(legacy_table[1]) / sum(legacy_table) else NA_real_
  marker_purity <- as.numeric(marker_table[1]) / sum(marker_table)
  final_label <- if (
    length(legacy_table) && sum(legacy_table) >= 50 && legacy_purity >= 0.55
  ) legacy_majority else marker_majority
  data.frame(
    cluster_rna = cluster_id,
    cells = length(index),
    cells_with_legacy_label = length(legacy_labels),
    legacy_majority = legacy_majority,
    legacy_purity = legacy_purity,
    marker_majority = marker_majority,
    marker_purity = marker_purity,
    final_broad_celltype = final_label,
    legacy_marker_majority_agree = !is.na(legacy_majority) && legacy_majority == marker_majority,
    stringsAsFactors = FALSE
  )
})
annotation_summary <- do.call(rbind, annotation_rows)
write.csv(
  annotation_summary,
  file.path(out_dir, "cluster_broad_annotation_summary.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

cluster_to_label <- setNames(
  annotation_summary$final_broad_celltype,
  annotation_summary$cluster_rna
)
object$celltype_broad_rechecked <- unname(cluster_to_label[as.character(object$cluster_rna)])

cluster_marker_scores <- aggregate(
  score_matrix,
  by = list(cluster_rna = object$cluster_rna),
  FUN = mean
)
write.csv(
  cluster_marker_scores,
  file.path(out_dir, "cluster_marker_module_scores.csv"),
  row.names = FALSE
)

celltype_sample_table <- as.data.frame.matrix(table(
  object$sample_id,
  object$celltype_broad_rechecked
))
celltype_sample_table$sample_id <- rownames(celltype_sample_table)
rownames(celltype_sample_table) <- NULL
celltype_sample_table <- celltype_sample_table[, c("sample_id", setdiff(colnames(celltype_sample_table), "sample_id"))]
write.csv(
  celltype_sample_table,
  file.path(out_dir, "broad_celltype_counts_by_sample.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

legacy_index <- !is.na(object$previous_celltype_main)
confusion <- as.data.frame.matrix(table(
  previous = object$previous_celltype_main[legacy_index],
  rechecked = object$celltype_broad_rechecked[legacy_index]
))
confusion$previous_celltype_main <- rownames(confusion)
rownames(confusion) <- NULL
confusion <- confusion[, c("previous_celltype_main", setdiff(colnames(confusion), "previous_celltype_main"))]
write.csv(
  confusion,
  file.path(out_dir, "previous_vs_rechecked_celltype_confusion.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

plot_celltype <- DimPlot(
  object,
  reduction = "umap.rpca",
  group.by = "celltype_broad_rechecked",
  raster = TRUE,
  label = TRUE,
  repel = TRUE
) + ggtitle("RNA-only rechecked broad cell types") + NoLegend()
plot_sample <- DimPlot(
  object,
  reduction = "umap.rpca",
  group.by = "sample_id",
  raster = TRUE
) + ggtitle("Samples after RPCA integration")
plot_severity <- DimPlot(
  object,
  reduction = "umap.rpca",
  group.by = "severity",
  raster = TRUE
) + ggtitle("Severity groups")

ggsave(
  file.path(out_dir, "umap_celltype_sample_severity.pdf"),
  plot_celltype + plot_sample + plot_severity,
  width = 18,
  height = 6,
  useDingbats = FALSE
)
ggsave(
  file.path(out_dir, "umap_celltype_sample_severity.png"),
  plot_celltype + plot_sample + plot_severity,
  width = 18,
  height = 6,
  dpi = 180
)

cell_meta <- object[[]]
cell_meta$cell_id <- rownames(cell_meta)
metadata_columns <- c(
  "cell_id", "raw_barcode", "sample_id", "severity", "allergy_group",
  "plastic_group", "batch_wave", "nCount_RNA", "nFeature_RNA",
  "percent.mt", "scrublet_score", "in_legacy_RDS", "previous_celltype_main",
  "previous_celltype_fine", "cluster_rna", "marker_broad_call",
  "celltype_broad_rechecked", "legacy_annotation_usable", "ADT_matrix_usable"
)
cell_con <- gzfile(file.path(out_dir, "cell_metadata_clustered.tsv.gz"), open = "wt")
write.table(
  cell_meta[, metadata_columns],
  file = cell_con,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE,
  na = "NA"
)
close(cell_con)

output_path <- file.path(out_dir, "23samples_RNA_reclustered_annotated.rds")
message("Saving reclustered RNA object: ", output_path)
saveRDS(object, output_path, compress = "gzip")

agreement_rate <- mean(
  object$previous_celltype_main[legacy_index] == object$celltype_broad_rechecked[legacy_index]
)
report_lines <- c(
  "# 23例 RNA 重聚类与宽谱细胞类型复核",
  "",
  paste0("- 运行时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("- QC后细胞数：", format(ncol(object), big.mark = ",")),
  paste0("- 使用可变基因数：", length(variable_features)),
  "- 标准化：LogNormalize。",
  "- 批次整合：Seurat RPCA；整合结果只用于邻居图、聚类和UMAP。",
  "- RPCA参考样本：V01、V18、V07、V14，覆盖两个采样波次及轻重症/过敏/塑形组合，减少23批次全配对整合的计算量。",
  "- 差异表达不会使用 RPCA 整合值，将使用 RNA 原始计数做样本级 pseudobulk。",
  paste0("- 0.6分辨率聚类数：", length(unique(object$cluster_rna))),
  paste0("- 与旧主类注释在可重叠细胞中的一致率：", sprintf("%.2f%%", 100 * agreement_rate)),
  "",
  "宽谱细胞类型以新聚类为单位确定：旧注释覆盖充分且纯度至少55%的聚类采用旧主类多数标签，否则采用独立标志基因模块得分多数标签。`cluster_broad_annotation_summary.csv` 保留了两套证据和一致性，后续亚群分析不会直接照搬旧细分标签。"
)
writeLines(report_lines, file.path(out_dir, "RNA重聚类与细胞类型复核报告.md"), useBytes = TRUE)

session_con <- file(file.path(out_dir, "sessionInfo.txt"), open = "wt")
writeLines(capture.output(sessionInfo()), session_con)
close(session_con)

message("RNA clustering and broad annotation completed")
