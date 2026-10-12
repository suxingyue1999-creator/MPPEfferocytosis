#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(edgeR)
})

set.seed(20260721)

workspace <- Sys.getenv("MPP_WORKSPACE", "inputs/workspace")
work_root <- file.path(workspace, "工作文件", "23例单细胞重新分析")
input_path <- file.path(
  work_root,
  "results",
  "02_rna_reanalysis",
  "23samples_RNA_reclustered_annotated.rds"
)
out_dir <- file.path(work_root, "results", "03_pseudobulk")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_path)) stop("Clustered RNA object is missing: ", input_path)
message("Loading clustered RNA object")
object <- readRDS(input_path)
DefaultAssay(object) <- "RNA"
if (!"celltype_broad_rechecked" %in% colnames(object[[]])) {
  stop("celltype_broad_rechecked is absent")
}
if (length(Layers(object[["RNA"]], search = "^counts$")) != 1L) {
  object <- JoinLayers(object, assay = "RNA")
}
counts <- LayerData(object[["RNA"]], layer = "counts")

cell_meta <- object[[]]
group_id <- paste(cell_meta$sample_id, cell_meta$celltype_broad_rechecked, sep = "||")
group_factor <- factor(group_id, levels = unique(group_id))
design_cells <- sparse.model.matrix(~ 0 + group_factor)
colnames(design_cells) <- levels(group_factor)
message("Aggregating raw counts to sample x broad-cell-type pseudobulks")
pseudobulk_counts <- counts %*% design_cells
rm(counts, design_cells)
gc(verbose = FALSE)

group_parts <- strsplit(colnames(pseudobulk_counts), "\\|\\|")
pseudobulk_meta <- data.frame(
  pseudobulk_id = colnames(pseudobulk_counts),
  sample_id = vapply(group_parts, `[[`, character(1), 1L),
  celltype = vapply(group_parts, `[[`, character(1), 2L),
  stringsAsFactors = FALSE
)
pseudobulk_meta$cell_count <- as.integer(table(group_factor)[pseudobulk_meta$pseudobulk_id])

sample_columns <- c(
  "sample_id", "severity", "allergy_group", "plastic_group", "batch_wave"
)
sample_meta <- unique(cell_meta[, sample_columns])
if (anyDuplicated(sample_meta$sample_id)) stop("Sample-level metadata is not unique")
pseudobulk_meta <- merge(
  pseudobulk_meta,
  sample_meta,
  by = "sample_id",
  all.x = TRUE,
  sort = FALSE
)
pseudobulk_meta <- pseudobulk_meta[match(colnames(pseudobulk_counts), pseudobulk_meta$pseudobulk_id), ]
if (anyNA(pseudobulk_meta$severity)) stop("Pseudobulk metadata mapping failed")

write.csv(
  pseudobulk_meta,
  file.path(out_dir, "pseudobulk_sample_metadata.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
saveRDS(
  pseudobulk_counts,
  file.path(out_dir, "pseudobulk_raw_counts.rds"),
  compress = "gzip"
)

contrast_map <- c(
  severity = "severity重症",
  allergy = "allergy_group是",
  plastic = "plastic_group是"
)
contrast_labels <- c(
  severity = "重症_vs_轻症",
  allergy = "过敏_vs_不过敏",
  plastic = "塑形_vs_非塑形"
)

result_summary <- list()
summary_index <- 0L
celltypes <- sort(unique(pseudobulk_meta$celltype))
for (celltype in celltypes) {
  message("Pseudobulk DE: ", celltype)
  keep_columns <- which(pseudobulk_meta$celltype == celltype & pseudobulk_meta$cell_count >= 20)
  current_meta <- pseudobulk_meta[keep_columns, , drop = FALSE]
  current_counts <- pseudobulk_counts[, keep_columns, drop = FALSE]

  current_meta$batch_wave <- factor(current_meta$batch_wave, levels = c("2023_Dec", "2024_Spring"))
  current_meta$severity <- factor(current_meta$severity, levels = c("轻症", "重症"))
  current_meta$allergy_group <- factor(current_meta$allergy_group, levels = c("否", "是"))
  current_meta$plastic_group <- factor(current_meta$plastic_group, levels = c("否", "是"))
  design <- model.matrix(
    ~ batch_wave + severity + allergy_group + plastic_group,
    data = current_meta
  )
  design_rank <- qr(design)$rank
  design_full_rank <- design_rank == ncol(design)

  if (nrow(current_meta) < 8L || !design_full_rank) {
    for (contrast_name in names(contrast_map)) {
      summary_index <- summary_index + 1L
      result_summary[[summary_index]] <- data.frame(
        celltype = celltype,
        contrast = contrast_labels[[contrast_name]],
        samples_used = nrow(current_meta),
        min_cells_per_pseudobulk = if (nrow(current_meta)) min(current_meta$cell_count) else NA,
        design_rank = design_rank,
        design_columns = ncol(design),
        genes_tested = NA_integer_,
        FDR_lt_0_05 = NA_integer_,
        FDR_lt_0_05_abs_logFC_ge_1 = NA_integer_,
        status = "skipped_insufficient_samples_or_rank",
        stringsAsFactors = FALSE
      )
    }
    next
  }

  dge <- DGEList(counts = current_counts)
  keep_genes <- filterByExpr(dge, design = design)
  dge <- dge[keep_genes, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge, method = "TMM")
  dge <- estimateDisp(dge, design = design, robust = TRUE)
  fit <- glmQLFit(dge, design = design, robust = TRUE)

  celltype_dir <- file.path(out_dir, gsub("[^A-Za-z0-9_]+", "_", celltype))
  dir.create(celltype_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(
    data.frame(sample = current_meta$sample_id, design, check.names = FALSE),
    file.path(celltype_dir, "design_matrix.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )

  for (contrast_name in names(contrast_map)) {
    coefficient <- contrast_map[[contrast_name]]
    label <- contrast_labels[[contrast_name]]
    if (!coefficient %in% colnames(design)) {
      summary_index <- summary_index + 1L
      result_summary[[summary_index]] <- data.frame(
        celltype = celltype,
        contrast = label,
        samples_used = nrow(current_meta),
        min_cells_per_pseudobulk = min(current_meta$cell_count),
        design_rank = design_rank,
        design_columns = ncol(design),
        genes_tested = nrow(dge),
        FDR_lt_0_05 = NA_integer_,
        FDR_lt_0_05_abs_logFC_ge_1 = NA_integer_,
        status = "skipped_coefficient_absent",
        stringsAsFactors = FALSE
      )
      next
    }
    test <- glmQLFTest(fit, coef = coefficient)
    result <- topTags(test, n = Inf, sort.by = "PValue")$table
    result$gene <- rownames(result)
    result <- result[, c("gene", setdiff(colnames(result), "gene"))]
    result$significant_FDR_0_05 <- result$FDR < 0.05
    result$significant_FDR_0_05_abs_logFC_1 <- result$FDR < 0.05 & abs(result$logFC) >= 1

    output_stem <- paste0(contrast_name, "_", gsub("[^A-Za-z0-9_]+", "_", label))
    write.csv(
      result,
      file.path(celltype_dir, paste0(output_stem, "_all_genes.csv")),
      row.names = FALSE,
      fileEncoding = "UTF-8"
    )
    write.csv(
      result[result$significant_FDR_0_05_abs_logFC_1, ],
      file.path(celltype_dir, paste0(output_stem, "_FDR005_logFC1.csv")),
      row.names = FALSE,
      fileEncoding = "UTF-8"
    )

    summary_index <- summary_index + 1L
    result_summary[[summary_index]] <- data.frame(
      celltype = celltype,
      contrast = label,
      samples_used = nrow(current_meta),
      min_cells_per_pseudobulk = min(current_meta$cell_count),
      design_rank = design_rank,
      design_columns = ncol(design),
      genes_tested = nrow(result),
      FDR_lt_0_05 = sum(result$significant_FDR_0_05),
      FDR_lt_0_05_abs_logFC_ge_1 = sum(result$significant_FDR_0_05_abs_logFC_1),
      status = "completed",
      stringsAsFactors = FALSE
    )
  }
}

result_summary <- do.call(rbind, result_summary)
write.csv(
  result_summary,
  file.path(out_dir, "pseudobulk_DE_summary.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

group_counts <- do.call(rbind, lapply(c("severity", "allergy_group", "plastic_group"), function(variable) {
  tab <- table(sample_meta$batch_wave, sample_meta[[variable]])
  out <- as.data.frame(tab, stringsAsFactors = FALSE)
  colnames(out) <- c("batch_wave", "group_level", "sample_count")
  out$group_variable <- variable
  out[, c("group_variable", "batch_wave", "group_level", "sample_count")]
}))
write.csv(
  group_counts,
  file.path(out_dir, "clinical_group_by_batch_wave.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

report_lines <- c(
  "# 23例样本级 pseudobulk 差异分析",
  "",
  paste0("- 运行时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("- 宽谱细胞类型数：", length(celltypes)),
  "- 计数输入：RNA 原始 UMI counts，不使用 RPCA 或其他校正表达矩阵。",
  "- 统计单位：患者样本；每个样本与细胞类型组合至少20个细胞才进入模型。",
  "- 模型：`~ batch_wave + severity + allergy_group + plastic_group`。",
  "- 批次变量：采样时间分为 2023年12月与 2024年春季两波。精确日期有12个水平且部分只有1例，直接纳入会消耗过多自由度，因此不用精确日期作固定效应。",
  "- 检验：edgeR quasi-likelihood，TMM标准化，robust dispersion。",
  "- 主要显著标准：FDR < 0.05 且 |log2FC| >= 1。",
  "",
  "`pseudobulk_DE_summary.csv` 给出每个细胞类型和对比的样本数、设计矩阵秩、检验基因数及显著基因数。"
)
writeLines(report_lines, file.path(out_dir, "pseudobulk分析说明.md"), useBytes = TRUE)

session_con <- file(file.path(out_dir, "sessionInfo.txt"), open = "wt")
writeLines(capture.output(sessionInfo()), session_con)
close(session_con)

message("Pseudobulk DE completed")
