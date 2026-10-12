#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
})

workspace <- Sys.getenv("MPP_WORKSPACE", "inputs/workspace")
input_root <- file.path(workspace, "scRNA-seq_data260717")
work_root <- file.path(workspace, "工作文件", "23例单细胞重新分析")
out_dir <- file.path(work_root, "results", "00_audit")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

object_path <- file.path(
  input_root,
  "WTD2024-04-006_23samples",
  "combin.data.RDS"
)
clinical_path <- file.path(workspace, "工作文件", "单细胞23例样本临床信息.csv")
rna_root <- file.path(
  input_root,
  "WTD2024-04-006(RNA&ADT report)",
  "WTD2024-04-006(RNA report)",
  "0_Rhapsody"
)
report_root <- dirname(rna_root)
scrublet_root <- file.path(
  report_root,
  "1_data_qc_statistics",
  "1_2_scrublet"
)
annotation_main_path <- file.path(
  report_root,
  "7_celltype_annotation",
  "7_1_SingleR",
  "cell_sample_cluster_annotation_RNA_type.main.xls"
)
annotation_fine_path <- file.path(
  report_root,
  "7_celltype_annotation",
  "7_1_SingleR",
  "cell_sample_cluster_annotation_RNA_type.fine.xls"
)

required_paths <- c(
  object_path,
  clinical_path,
  rna_root,
  scrublet_root,
  annotation_main_path,
  annotation_fine_path
)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths)) {
  stop("Missing required input(s):\n", paste(missing_paths, collapse = "\n"))
}

read_utf8_csv <- function(path) {
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
}

read_gzip_lines <- function(path, n = -1L) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  readLines(con, n = n, warn = FALSE)
}

read_mtx_dimensions <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(line)) stop("No Matrix Market dimensions found: ", path)
    if (!startsWith(line, "%")) break
  }
  values <- as.numeric(strsplit(trimws(line), "[[:space:]]+")[[1]])
  if (length(values) != 3L) stop("Invalid Matrix Market dimensions: ", path)
  setNames(values, c("features", "barcodes", "nonzero_entries"))
}

safe_median <- function(x) median(x, na.rm = TRUE)
safe_mad <- function(x) mad(x, constant = 1.4826, na.rm = TRUE)

message("Loading legacy Seurat object (read-only): ", object_path)
legacy <- readRDS(object_path)
if (!inherits(legacy, "Seurat")) stop("RDS is not a Seurat object")
if (!"RNA" %in% names(legacy@assays)) stop("RNA assay is absent")

counts <- legacy@assays$RNA@counts
meta <- legacy@meta.data
if (!identical(colnames(counts), rownames(meta))) {
  stop("Counts columns and metadata row names are not identical")
}

clinical <- read_utf8_csv(clinical_path)
if (anyDuplicated(clinical$sample_code)) stop("Clinical sample_code contains duplicates")

object_samples <- sort(unique(meta$orig.ident))
clinical_samples <- sort(unique(clinical$sample_code))
raw_samples <- sort(list.dirs(rna_root, recursive = FALSE, full.names = FALSE))
expected_samples <- sort(intersect(raw_samples, clinical_samples))

sample_set_audit <- data.frame(
  sample_id = sort(unique(c(object_samples, clinical_samples, raw_samples))),
  stringsAsFactors = FALSE
)
sample_set_audit$in_object <- sample_set_audit$sample_id %in% object_samples
sample_set_audit$in_clinical <- sample_set_audit$sample_id %in% clinical_samples
sample_set_audit$in_raw_matrix <- sample_set_audit$sample_id %in% raw_samples
write.csv(
  sample_set_audit,
  file.path(out_dir, "sample_set_audit.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
if (!all(sample_set_audit$in_object & sample_set_audit$in_clinical & sample_set_audit$in_raw_matrix)) {
  stop("Sample sets differ; inspect sample_set_audit.csv")
}

cell_meta <- meta
cell_meta$cell_id <- rownames(cell_meta)
cell_meta$sample_id <- cell_meta$orig.ident
cell_meta$raw_barcode <- sub("_[0-9]+$", "", cell_meta$cell_id)
cell_meta$mito_percent_value <- if (max(cell_meta$mito.percent, na.rm = TRUE) <= 1) {
  100 * cell_meta$mito.percent
} else {
  cell_meta$mito.percent
}

clinical_columns <- c(
  "sample_code", "severity",
  "allergy_group", "plastic_group", "age", "sex"
)
clinical_index <- match(cell_meta$sample_id, clinical$sample_code)
if (anyNA(clinical_index)) stop("Some object cells do not map to clinical samples")
for (column in setdiff(clinical_columns, "sample_code")) {
  cell_meta[[column]] <- clinical[[column]][clinical_index]
}
cell_meta$legacy_annotation_usable <- cell_meta$sample_id != "V20"
cell_meta$ADT_matrix_usable <- TRUE
cell_meta$RNA_analysis_usable <- TRUE

main_annotation <- read.delim(
  annotation_main_path,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "\""
)
fine_annotation <- read.delim(
  annotation_fine_path,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "\""
)
main_index <- match(cell_meta$cell_id, main_annotation$barcodes)
fine_index <- match(cell_meta$cell_id, fine_annotation$barcodes)
cell_meta$previous_celltype_main <- main_annotation$cluster.combin.labels[main_index]
cell_meta$previous_celltype_fine <- fine_annotation$cluster.combin.labels[fine_index]

cell_meta$scrublet_score <- NA_real_
cell_meta$scrublet_prediction <- NA
scrublet_raw_rows <- vector("list", length(object_samples))
names(scrublet_raw_rows) <- object_samples
for (sample_id in object_samples) {
  path <- file.path(scrublet_root, paste0(sample_id, ".doublet_prediction.txt"))
  if (!file.exists(path)) stop("Missing Scrublet file: ", path)
  scrub <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  scrub$prediction <- tolower(as.character(scrub$prediction)) == "true"
  cell_index <- which(cell_meta$sample_id == sample_id)
  scrub_index <- match(cell_meta$raw_barcode[cell_index], scrub$barcodes)
  cell_meta$scrublet_score[cell_index] <- scrub$score[scrub_index]
  cell_meta$scrublet_prediction[cell_index] <- scrub$prediction[scrub_index]
  scrublet_raw_rows[[sample_id]] <- data.frame(
    sample_id = sample_id,
    raw_scrublet_cells = nrow(scrub),
    raw_scrublet_doublets = sum(scrub$prediction, na.rm = TRUE),
    raw_scrublet_doublet_rate = mean(scrub$prediction, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}
scrublet_raw_summary <- do.call(rbind, scrublet_raw_rows)

qc_by_sample <- split(seq_len(nrow(cell_meta)), cell_meta$sample_id)
threshold_rows <- lapply(names(qc_by_sample), function(sample_id) {
  index <- qc_by_sample[[sample_id]]
  nf <- cell_meta$nFeature_RNA[index]
  nc <- cell_meta$nCount_RNA[index]
  mt <- cell_meta$mito_percent_value[index]
  data.frame(
    sample_id = sample_id,
    min_features = max(200, floor(safe_median(nf) - 3 * safe_mad(nf))),
    max_features = min(8000, ceiling(safe_median(nf) + 3 * safe_mad(nf))),
    min_counts = max(500, floor(safe_median(nc) - 3 * safe_mad(nc))),
    max_counts = ceiling(safe_median(nc) + 3 * safe_mad(nc)),
    max_mito_percent = min(25, safe_median(mt) + 3 * safe_mad(mt)),
    stringsAsFactors = FALSE
  )
})
thresholds <- do.call(rbind, threshold_rows)
write.csv(
  thresholds,
  file.path(out_dir, "provisional_qc_thresholds_by_sample.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

threshold_index <- match(cell_meta$sample_id, thresholds$sample_id)
cell_meta$qc_gene_count_pass <-
  cell_meta$nFeature_RNA >= thresholds$min_features[threshold_index] &
  cell_meta$nFeature_RNA <= thresholds$max_features[threshold_index]
cell_meta$qc_umi_count_pass <-
  cell_meta$nCount_RNA >= thresholds$min_counts[threshold_index] &
  cell_meta$nCount_RNA <= thresholds$max_counts[threshold_index]
cell_meta$qc_mito_pass <-
  cell_meta$mito_percent_value <= thresholds$max_mito_percent[threshold_index]
cell_meta$qc_doublet_pass <- !is.na(cell_meta$scrublet_prediction) & !cell_meta$scrublet_prediction
cell_meta$provisional_qc_pass <- with(
  cell_meta,
  qc_gene_count_pass & qc_umi_count_pass & qc_mito_pass & qc_doublet_pass
)

raw_audit_rows <- lapply(expected_samples, function(sample_id) {
  sample_dir <- file.path(rna_root, sample_id)
  barcode_path <- file.path(sample_dir, "barcodes.tsv.gz")
  feature_path <- file.path(sample_dir, "features.tsv.gz")
  matrix_path <- file.path(sample_dir, "matrix.mtx.gz")
  dimensions <- read_mtx_dimensions(matrix_path)
  raw_barcodes <- read_gzip_lines(barcode_path)
  raw_features <- read_gzip_lines(feature_path)
  object_barcodes <- cell_meta$raw_barcode[cell_meta$sample_id == sample_id]
  overlap_count <- sum(object_barcodes %in% raw_barcodes)
  data.frame(
    sample_id = sample_id,
    raw_barcode_lines = length(raw_barcodes),
    raw_feature_lines = length(raw_features),
    matrix_features = dimensions[["features"]],
    matrix_barcodes = dimensions[["barcodes"]],
    matrix_nonzero_entries = dimensions[["nonzero_entries"]],
    object_cells = length(object_barcodes),
    overlapping_object_cells = overlap_count,
    object_cells_are_raw_subset = all(object_barcodes %in% raw_barcodes),
    object_overlap_rate = overlap_count / length(object_barcodes),
    barcodes_identical_as_set = setequal(raw_barcodes, object_barcodes),
    matrix_size_bytes = file.info(matrix_path)$size,
    stringsAsFactors = FALSE
  )
})
raw_audit <- do.call(rbind, raw_audit_rows)
write.csv(
  raw_audit,
  file.path(out_dir, "raw_matrix_audit.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

sample_summary_rows <- lapply(names(qc_by_sample), function(sample_id) {
  index <- qc_by_sample[[sample_id]]
  cindex <- match(sample_id, clinical$sample_code)
  sindex <- match(sample_id, scrublet_raw_summary$sample_id)
  data.frame(
    sample_id = sample_id,
    severity = clinical$severity[cindex],
    allergy_group = clinical$allergy_group[cindex],
    plastic_group = clinical$plastic_group[cindex],
    cells_before_qc = length(index),
    median_features = safe_median(cell_meta$nFeature_RNA[index]),
    median_umi = safe_median(cell_meta$nCount_RNA[index]),
    median_mito_percent = safe_median(cell_meta$mito_percent_value[index]),
    raw_scrublet_cells = scrublet_raw_summary$raw_scrublet_cells[sindex],
    raw_scrublet_doublets = scrublet_raw_summary$raw_scrublet_doublets[sindex],
    raw_scrublet_doublet_rate = scrublet_raw_summary$raw_scrublet_doublet_rate[sindex],
    legacy_object_scrublet_doublets = sum(cell_meta$scrublet_prediction[index], na.rm = TRUE),
    legacy_object_scrublet_unmapped = sum(is.na(cell_meta$scrublet_prediction[index])),
    provisional_qc_pass = sum(cell_meta$provisional_qc_pass[index], na.rm = TRUE),
    provisional_qc_pass_rate = mean(cell_meta$provisional_qc_pass[index], na.rm = TRUE),
    previous_main_annotation_missing = sum(is.na(cell_meta$previous_celltype_main[index])),
    previous_fine_annotation_missing = sum(is.na(cell_meta$previous_celltype_fine[index])),
    legacy_annotation_usable = sample_id != "V20",
    ADT_matrix_usable = TRUE,
    stringsAsFactors = FALSE
  )
})
sample_summary <- do.call(rbind, sample_summary_rows)
sample_summary <- sample_summary[match(clinical$sample_code, sample_summary$sample_id), ]
write.csv(
  sample_summary,
  file.path(out_dir, "sample_qc_audit.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

cell_output_columns <- c(
  "cell_id", "raw_barcode", "sample_id", "severity", "allergy_group", "plastic_group", "age", "sex", "nCount_RNA", "nFeature_RNA",
  "mito_percent_value", "ribo.percent", "diss.percent", "seurat_clusters",
  "previous_celltype_main", "previous_celltype_fine", "scrublet_score",
  "scrublet_prediction", "qc_gene_count_pass", "qc_umi_count_pass",
  "qc_mito_pass", "qc_doublet_pass", "provisional_qc_pass",
  "RNA_analysis_usable", "legacy_annotation_usable", "ADT_matrix_usable"
)
cell_output <- cell_meta[, cell_output_columns]
cell_con <- gzfile(file.path(out_dir, "cell_metadata_audit.tsv.gz"), open = "wt")
write.table(
  cell_output,
  file = cell_con,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE,
  na = "NA"
)
close(cell_con)

group_summary <- do.call(rbind, lapply(c("severity", "allergy_group", "plastic_group"), function(group_var) {
  tab <- table(clinical[[group_var]], useNA = "ifany")
  data.frame(
    group_variable = group_var,
    group_level = names(tab),
    sample_count = as.integer(tab),
    stringsAsFactors = FALSE
  )
}))
write.csv(
  group_summary,
  file.path(out_dir, "clinical_group_counts.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

object_version <- tryCatch(as.character(legacy@version), error = function(e) "unknown")
audit_lines <- c(
  "# 23例 RNA 数据首轮审计",
  "",
  paste0("- 审计时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("- 原始 RDS：`", object_path, "`（只读）"),
  paste0("- 对象保存时 Seurat 版本：", object_version),
  paste0("- 当前 R：", R.version.string),
  paste0("- 当前 Seurat：", if (requireNamespace("Seurat", quietly = TRUE)) as.character(packageVersion("Seurat")) else "未安装"),
  paste0("- RNA 基因数：", format(nrow(counts), big.mark = ",")),
  paste0("- 细胞数：", format(ncol(counts), big.mark = ",")),
  paste0("- 样本数：", length(object_samples)),
  paste0("- 原始计数保留：", inherits(counts, "sparseMatrix")),
  paste0("- 对象/临床/原始矩阵样本集合完全一致：", all(sample_set_audit$in_object & sample_set_audit$in_clinical & sample_set_audit$in_raw_matrix)),
  paste0("- 旧 RDS 细胞是当前矩阵条码严格子集的样本数：", sum(raw_audit$object_cells_are_raw_subset), "/23"),
  paste0("- V20 旧 RDS/当前矩阵重叠细胞数：", raw_audit$overlapping_object_cells[raw_audit$sample_id == "V20"], "/", raw_audit$object_cells[raw_audit$sample_id == "V20"]),
  paste0("- Scrublet 未映射细胞数：", sum(is.na(cell_meta$scrublet_prediction))),
  paste0("- 原主类注释缺失细胞数：", sum(is.na(cell_meta$previous_celltype_main))),
  paste0("- 原细分注释缺失细胞数：", sum(is.na(cell_meta$previous_celltype_fine))),
  "",
  "## 已固定的分析规则",
  "",
  "- Public inputs use coded sample identifiers consistently.",
  "- 轻症/重症、过敏/不过敏、塑形/非塑形均来自顶层临床分组表。",
  "- 23例和10例是独立队列，不做患者级配对。",
  "- 当前23个 `0_Rhapsody` 矩阵均在同一文件内包含 RNA 和 ADT；`ADT_matrix_usable=TRUE`。历史字段 `legacy_annotation_usable=FALSE` 仅表示旧供应商 RDS/WNN 标签不可直接复用。",
  "- 旧 RDS 的 V20 与当前矩阵不是同一细胞版本。正式 RNA/ADT 重分析从当前23个 `0_Rhapsody` 矩阵起步；旧 RDS 只用于历史结果对照。",
  "- 批次校正仅用于可视化和聚类；pseudobulk 使用原始计数，并在设计矩阵中处理可识别的批次或日期因素。",
  "",
  "## QC 说明",
  "",
  "`provisional_qc_thresholds_by_sample.csv` 是基于旧 RDS 中每个样本中位数和 MAD 的参考阈值，仅用于审计。正式过滤从当前原始矩阵重新计算阈值，尤其要单独检查 V20，避免用全队列统一阈值误删。",
  "",
  "原始 Scrublet 文件中存在预测双细胞；旧 RDS 内为0，是因为供应商已在合并对象前将预测双细胞过滤掉。正式重分析会在当前原始矩阵层面重新合并这些 Scrublet 标记。",
  "",
  "## 输出",
  "",
  "- `sample_set_audit.csv`：三处样本集合核对。",
  "- `raw_matrix_audit.csv`：每例矩阵维度、条码及对象一致性。",
  "- `sample_qc_audit.csv`：逐样本 QC、Scrublet 和首轮保留率。",
  "- `cell_metadata_audit.tsv.gz`：统一到细胞级的分组、旧注释、QC和双细胞信息。",
  "- `clinical_group_counts.csv`：三类临床分组样本数。"
)
writeLines(audit_lines, file.path(out_dir, "审计报告.md"), useBytes = TRUE)

session_con <- file(file.path(out_dir, "sessionInfo.txt"), open = "wt")
writeLines(capture.output(sessionInfo()), session_con)
close(session_con)

message("Audit completed: ", out_dir)
