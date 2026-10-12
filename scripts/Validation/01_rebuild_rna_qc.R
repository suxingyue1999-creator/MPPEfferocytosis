#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
})

set.seed(20260721)
options(future.globals.maxSize = 40 * 1024^3)

workspace <- Sys.getenv("MPP_WORKSPACE", "inputs/workspace")
input_root <- file.path(workspace, "scRNA-seq_data260717")
work_root <- file.path(workspace, "工作文件", "23例单细胞重新分析")
audit_dir <- file.path(work_root, "results", "00_audit")
out_dir <- file.path(work_root, "results", "01_qc")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

clinical_path <- file.path(workspace, "工作文件", "单细胞23例样本临床信息.csv")
rna_root <- file.path(
  input_root,
  "WTD2024-04-006(RNA&ADT report)",
  "WTD2024-04-006(RNA report)",
  "0_Rhapsody"
)
scrublet_root <- file.path(
  dirname(rna_root),
  "1_data_qc_statistics",
  "1_2_scrublet"
)
legacy_cell_meta_path <- file.path(audit_dir, "cell_metadata_audit.tsv.gz")

required_paths <- c(clinical_path, rna_root, scrublet_root, legacy_cell_meta_path)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths)) {
  stop("Missing required input(s):\n", paste(missing_paths, collapse = "\n"))
}

clinical <- read.csv(
  clinical_path,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  fileEncoding = "UTF-8-BOM"
)
samples <- clinical$sample_code
if (length(samples) != 23L || anyDuplicated(samples)) {
  stop("Expected 23 unique clinical sample codes")
}

legacy_cell_meta <- read.delim(
  gzfile(legacy_cell_meta_path),
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = "NA"
)
legacy_key <- paste(legacy_cell_meta$sample_id, legacy_cell_meta$raw_barcode, sep = "::")
if (anyDuplicated(legacy_key)) stop("Legacy sample/barcode keys are not unique")

safe_median <- function(x) median(x, na.rm = TRUE)
safe_mad <- function(x) mad(x, constant = 1.4826, na.rm = TRUE)

sample_objects <- vector("list", length(samples))
names(sample_objects) <- samples
qc_rows <- vector("list", length(samples))
names(qc_rows) <- samples
threshold_rows <- vector("list", length(samples))
names(threshold_rows) <- samples

for (sample_id in samples) {
  message("Reading current RNA matrix: ", sample_id)
  sample_dir <- file.path(rna_root, sample_id)
  matrix_paths <- file.path(sample_dir, c("matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz"))
  if (!all(file.exists(matrix_paths))) {
    stop("Incomplete matrix directory: ", sample_dir)
  }

  counts <- Read10X(
    data.dir = sample_dir,
    gene.column = 2,
    unique.features = TRUE,
    strip.suffix = FALSE
  )
  if (is.list(counts)) {
    if (!"Gene Expression" %in% names(counts)) {
      stop("Gene Expression matrix absent for ", sample_id)
    }
    counts <- counts[["Gene Expression"]]
  }
  object <- CreateSeuratObject(
    counts = counts,
    project = sample_id,
    assay = "RNA",
    min.cells = 0,
    min.features = 0
  )
  rm(counts)
  object <- RenameCells(object, add.cell.id = sample_id)
  object$sample_id <- sample_id
  object$raw_barcode <- sub(paste0("^", sample_id, "_"), "", Cells(object))

  clinical_index <- match(sample_id, clinical$sample_code)
  object$severity <- clinical$severity[clinical_index]
  object$allergy_group <- clinical$allergy_group[clinical_index]
  object$plastic_group <- clinical$plastic_group[clinical_index]
  object$batch_wave <- as.character(clinical$batch_wave[clinical_index])
  object$age <- clinical$age[clinical_index]
  object$sex <- clinical$sex[clinical_index]
  object$RNA_analysis_usable <- TRUE
  object$legacy_annotation_usable <- sample_id != "V20"
  object$ADT_matrix_usable <- TRUE

  object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT-")
  object[["percent.ribo"]] <- PercentageFeatureSet(object, pattern = "^RP[SL]")

  scrublet_path <- file.path(scrublet_root, paste0(sample_id, ".doublet_prediction.txt"))
  scrublet <- read.csv(scrublet_path, stringsAsFactors = FALSE, check.names = FALSE)
  scrublet$prediction <- tolower(as.character(scrublet$prediction)) == "true"
  scrublet_index <- match(object$raw_barcode, scrublet$barcodes)
  if (anyNA(scrublet_index)) {
    stop(sample_id, ": current RNA barcodes do not all map to Scrublet results")
  }
  object$scrublet_score <- scrublet$score[scrublet_index]
  object$scrublet_prediction <- scrublet$prediction[scrublet_index]

  current_key <- paste(sample_id, object$raw_barcode, sep = "::")
  legacy_index <- match(current_key, legacy_key)
  object$in_legacy_RDS <- !is.na(legacy_index)
  object$previous_celltype_main <- legacy_cell_meta$previous_celltype_main[legacy_index]
  object$previous_celltype_fine <- legacy_cell_meta$previous_celltype_fine[legacy_index]
  object$previous_seurat_cluster <- legacy_cell_meta$seurat_clusters[legacy_index]

  metadata <- object[[]]
  median_features <- safe_median(metadata$nFeature_RNA)
  median_counts <- safe_median(metadata$nCount_RNA)
  median_mito <- safe_median(metadata$percent.mt)
  mad_features <- safe_mad(metadata$nFeature_RNA)
  mad_counts <- safe_mad(metadata$nCount_RNA)
  mad_mito <- safe_mad(metadata$percent.mt)

  min_features <- 200
  max_features <- min(8000, ceiling(median_features + 5 * mad_features))
  min_counts <- 500
  max_counts <- ceiling(median_counts + 5 * mad_counts)
  max_mito <- min(25, median_mito + 5 * mad_mito)
  if (!is.finite(max_features) || max_features <= min_features) max_features <- 8000
  if (!is.finite(max_counts) || max_counts <= min_counts) max_counts <- max(metadata$nCount_RNA)
  if (!is.finite(max_mito) || max_mito <= 0) max_mito <- 25

  object$qc_min_features_pass <- object$nFeature_RNA >= min_features
  object$qc_max_features_pass <- object$nFeature_RNA <= max_features
  object$qc_min_counts_pass <- object$nCount_RNA >= min_counts
  object$qc_max_counts_pass <- object$nCount_RNA <= max_counts
  object$qc_mito_pass <- object$percent.mt <= max_mito
  object$qc_doublet_pass <- !object$scrublet_prediction
  object$qc_pass <- with(
    object[[]],
    qc_min_features_pass & qc_max_features_pass & qc_min_counts_pass &
      qc_max_counts_pass & qc_mito_pass & qc_doublet_pass
  )

  threshold_rows[[sample_id]] <- data.frame(
    sample_id = sample_id,
    min_features = min_features,
    max_features = max_features,
    min_counts = min_counts,
    max_counts = max_counts,
    max_mito_percent = max_mito,
    median_features_raw = median_features,
    median_counts_raw = median_counts,
    median_mito_percent_raw = median_mito,
    stringsAsFactors = FALSE
  )

  qc_rows[[sample_id]] <- data.frame(
    sample_id = sample_id,
    severity = object$severity[[1]],
    allergy_group = object$allergy_group[[1]],
    plastic_group = object$plastic_group[[1]],
    batch_wave = object$batch_wave[[1]],
    cells_raw_matrix = ncol(object),
    scrublet_doublets = sum(object$scrublet_prediction),
    cells_in_legacy_RDS = sum(object$in_legacy_RDS),
    fail_min_features = sum(!object$qc_min_features_pass),
    fail_max_features = sum(!object$qc_max_features_pass),
    fail_min_counts = sum(!object$qc_min_counts_pass),
    fail_max_counts = sum(!object$qc_max_counts_pass),
    fail_mito = sum(!object$qc_mito_pass),
    cells_qc_pass = sum(object$qc_pass),
    qc_pass_rate = mean(object$qc_pass),
    stringsAsFactors = FALSE
  )

  keep_cells <- Cells(object)[object$qc_pass]
  sample_objects[[sample_id]] <- subset(object, cells = keep_cells)
  rm(object, scrublet, metadata)
  gc(verbose = FALSE)
}

thresholds <- do.call(rbind, threshold_rows)
qc_summary <- do.call(rbind, qc_rows)
row.names(thresholds) <- NULL
row.names(qc_summary) <- NULL
write.csv(
  thresholds,
  file.path(out_dir, "qc_thresholds_by_sample.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  qc_summary,
  file.path(out_dir, "qc_summary_by_sample.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

message("Merging 23 QC-filtered sample objects")
merged <- merge(
  x = sample_objects[[1]],
  y = sample_objects[-1],
  merge.data = FALSE,
  merge.dr = FALSE
)
rm(sample_objects)
gc(verbose = FALSE)

if (anyDuplicated(Cells(merged))) stop("Merged object contains duplicated cell IDs")
if (ncol(merged) != sum(qc_summary$cells_qc_pass)) {
  stop("Merged cell count does not equal the sum of per-sample QC counts")
}

cell_meta <- merged[[]]
cell_meta$cell_id <- rownames(cell_meta)
cell_meta <- cell_meta[, c(
  "cell_id", "raw_barcode", "sample_id", "severity", "allergy_group", "plastic_group", "batch_wave",
  "age", "sex", "nCount_RNA", "nFeature_RNA",
  "percent.mt", "percent.ribo", "scrublet_score", "scrublet_prediction",
  "in_legacy_RDS", "previous_celltype_main", "previous_celltype_fine",
  "previous_seurat_cluster", "RNA_analysis_usable", "legacy_annotation_usable",
  "ADT_matrix_usable"
)]
cell_con <- gzfile(file.path(out_dir, "cell_metadata_after_qc.tsv.gz"), open = "wt")
write.table(
  cell_meta,
  file = cell_con,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE,
  na = "NA"
)
close(cell_con)

qc_long <- rbind(
  data.frame(sample_id = qc_summary$sample_id, metric = "Raw matrix", cells = qc_summary$cells_raw_matrix),
  data.frame(sample_id = qc_summary$sample_id, metric = "After QC", cells = qc_summary$cells_qc_pass)
)
qc_long$sample_id <- factor(qc_long$sample_id, levels = samples)
qc_long$metric <- factor(qc_long$metric, levels = c("Raw matrix", "After QC"))
plot <- ggplot(qc_long, aes(x = sample_id, y = cells, fill = metric)) +
  geom_col(position = "dodge", width = 0.8) +
  scale_fill_manual(values = c("Raw matrix" = "#6B7280", "After QC" = "#1B9E77")) +
  labs(x = NULL, y = "Cell count", fill = NULL, title = "Cells before and after RNA QC") +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "top")
ggsave(
  file.path(out_dir, "cells_before_after_qc.pdf"),
  plot,
  width = 10,
  height = 5.5,
  useDingbats = FALSE
)
ggsave(
  file.path(out_dir, "cells_before_after_qc.png"),
  plot,
  width = 10,
  height = 5.5,
  dpi = 180
)

object_path <- file.path(out_dir, "23samples_current_RNA_qc_counts.rds")
message("Saving QC counts object: ", object_path)
saveRDS(merged, object_path, compress = "gzip")

report_lines <- c(
  "# 23例当前 RNA 矩阵 QC 结果",
  "",
  paste0("- 运行时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("- 当前矩阵原始细胞数：", format(sum(qc_summary$cells_raw_matrix), big.mark = ",")),
  paste0("- Scrublet 预测双细胞数：", format(sum(qc_summary$scrublet_doublets), big.mark = ",")),
  paste0("- QC 后细胞数：", format(sum(qc_summary$cells_qc_pass), big.mark = ",")),
  paste0("- 总保留率：", sprintf("%.2f%%", 100 * sum(qc_summary$cells_qc_pass) / sum(qc_summary$cells_raw_matrix))),
  paste0("- V20 当前矩阵细胞数：", qc_summary$cells_raw_matrix[qc_summary$sample_id == "V20"]),
  paste0("- V20 与旧 RDS 重叠细胞数：", qc_summary$cells_in_legacy_RDS[qc_summary$sample_id == "V20"]),
  "",
  "## 过滤规则",
  "",
  "- 每个样本至少200个基因、至少500 UMI。",
  "- 每个样本的高基因数和高 UMI 上限为中位数加5倍 MAD，并设置8,000基因硬上限。",
  "- 线粒体比例上限为样本中位数加5倍 MAD，最高不超过25%。",
  "- 去除供应商原始 Scrublet 文件中预测为双细胞的条码。",
  "- 阈值按样本单独计算，避免 V20 的低深度分布被其他样本阈值误伤。",
  "",
  "## 版本规则",
  "",
  "正式 RNA 重分析从当前23个 `0_Rhapsody` 矩阵起步。当前矩阵中的 RNA/ADT 同条码，`ADT_matrix_usable=TRUE`；旧 Seurat 3 RDS 仅用于历史注释和条码重叠对照，旧供应商 V20 WNN 标签不直接复用。"
)
writeLines(report_lines, file.path(out_dir, "RNA_QC报告.md"), useBytes = TRUE)

session_con <- file(file.path(out_dir, "sessionInfo.txt"), open = "wt")
writeLines(capture.output(sessionInfo()), session_con)
close(session_con)

message("RNA QC completed: ", out_dir)
