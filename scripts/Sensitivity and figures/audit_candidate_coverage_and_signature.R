#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(edgeR)
})

base <- Sys.getenv("MPP_WORKSPACE", "inputs/workspace")
article1 <- file.path(base, "工作文件", "文章1")
out_dir <- file.path(article1, "小验证_20260831")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

axis_path <- file.path(article1, "补充表", "S5_HIF1A候选节点.csv")
extra_path <- file.path(article1, "补充表", "S6_胞葬处理与旁路候选节点.csv")
meta_path <- file.path(
  base, "工作文件", "23例单细胞重新分析", "results", "07_lineage_pseudobulk",
  "Macrophage_Monocyte", "pseudobulk_sample_metadata.csv"
)
counts_path <- file.path(
  base, "工作文件", "23例单细胞重新分析", "results", "07_lineage_pseudobulk",
  "Macrophage_Monocyte", "pseudobulk_raw_counts.rds"
)
subtype_root <- file.path(
  base, "工作文件", "23例单细胞重新分析", "results", "07_lineage_pseudobulk",
  "Macrophage_Monocyte"
)

required <- c(axis_path, extra_path, meta_path, counts_path)
if (!all(file.exists(required))) {
  stop("Missing required input(s): ", paste(required[!file.exists(required)], collapse = "; "))
}

axis <- read.csv(axis_path, check.names = FALSE, fileEncoding = "UTF-8-BOM")
extra <- read.csv(extra_path, check.names = FALSE, fileEncoding = "UTF-8-BOM")
axis_genes <- unique(as.character(axis$gene))
extra_genes <- unique(as.character(extra$gene))
core_genes <- c("HIF1A", "THBS1", "CD36", "MERTK", "TGFB1", "MAPK1", "MAPK3", "SMAD2", "SMAD4")
panel <- unique(c(core_genes, axis_genes, extra_genes))
panel <- panel[!is.na(panel) & nzchar(panel)]

discovery <- axis[axis$gene %in% core_genes, c("gene", "log2FC", "PValue", "FDR"), drop = FALSE]
names(discovery) <- c("gene", "discovery_logFC", "discovery_PValue", "discovery_FDR")
discovery <- discovery[!duplicated(discovery$gene), , drop = FALSE]
discovery$discovery_direction <- ifelse(discovery$discovery_logFC > 0, "positive", ifelse(discovery$discovery_logFC < 0, "negative", "zero"))

meta <- read.csv(meta_path, check.names = FALSE, fileEncoding = "UTF-8-BOM")
counts <- readRDS(counts_path)
counts <- as.matrix(counts)
if (is.null(colnames(counts)) || is.null(rownames(counts))) stop("Pseudobulk counts lack row or column names")
meta <- meta[match(colnames(counts), meta$pseudobulk_id), , drop = FALSE]
if (anyNA(meta$pseudobulk_id)) stop("Pseudobulk metadata does not match count columns")

subtypes <- sort(unique(as.character(meta$subtype)))
rows <- list()
idx <- 0L

for (subtype in subtypes) {
  current_idx <- which(as.character(meta$subtype) == subtype)
  current <- meta[current_idx, , drop = FALSE]
  current_counts <- counts[, current_idx, drop = FALSE]
  eligible <- current$cell_count >= 20
  current <- current[eligible, , drop = FALSE]
  current_counts <- current_counts[, eligible, drop = FALSE]
  if (!nrow(current)) next

  current$batch_wave <- factor(current$batch_wave, levels = c("2023_Dec", "2024_Spring"))
  current$severity <- factor(current$severity, levels = c("轻症", "重症"))
  current$allergy_group <- factor(current$allergy_group, levels = c("否", "是"))
  current$plastic_group <- factor(current$plastic_group, levels = c("否", "是"))
  design <- model.matrix(~ batch_wave + severity + allergy_group + plastic_group, data = current)
  model_status <- if (nrow(current) < 8L || qr(design)$rank != ncol(design)) {
    "not_estimable_under_current_design"
  } else {
    "reported_from_existing_edgeR_output"
  }

  result_path <- file.path(subtype_root, subtype, "severity___vs___all_genes.csv")
  result <- if (file.exists(result_path)) read.csv(result_path, check.names = FALSE) else NULL
  mild_n <- sum(current$severity == "轻症")
  severe_n <- sum(current$severity == "重症")
  for (gene in panel) {
    idx <- idx + 1L
    hit <- if (!is.null(result)) result[result$gene == gene, , drop = FALSE] else result
    if (!is.null(hit) && nrow(hit) > 1L) hit <- hit[1, , drop = FALSE]
    discovery_row <- discovery[discovery$gene == gene, , drop = FALSE]
    feature_present <- gene %in% rownames(current_counts)
    expression_samples <- if (feature_present) sum(current_counts[gene, , drop = TRUE] > 0) else NA_integer_
    expansion_logfc <- if (!is.null(hit) && nrow(hit)) hit$logFC[[1]] else NA_real_
    expansion_p <- if (!is.null(hit) && nrow(hit)) hit$PValue[[1]] else NA_real_
    expansion_fdr <- if (!is.null(hit) && nrow(hit)) hit$FDR[[1]] else NA_real_
    expansion_direction <- if (is.na(expansion_logfc)) NA_character_ else if (expansion_logfc > 0) "positive" else if (expansion_logfc < 0) "negative" else "zero"
    discovery_direction <- if (nrow(discovery_row)) discovery_row$discovery_direction[[1]] else NA_character_
    direction_concordant <- if (is.na(discovery_direction) || is.na(expansion_direction)) NA else discovery_direction == expansion_direction
    rows[[idx]] <- data.frame(
      cohort = "23_patient_expanded",
      lineage = "Macrophage_Monocyte",
      subtype = subtype,
      samples_used = nrow(current),
      mild_n = mild_n,
      severe_n = severe_n,
      min_cells_per_pseudobulk = if (nrow(current)) min(current$cell_count) else NA_integer_,
      model_status = model_status,
      gene = gene,
      feature_present_in_count_matrix = feature_present,
      samples_with_positive_pseudobulk_count = expression_samples,
      discovery_logFC = if (nrow(discovery_row)) discovery_row$discovery_logFC[[1]] else NA_real_,
      discovery_PValue = if (nrow(discovery_row)) discovery_row$discovery_PValue[[1]] else NA_real_,
      discovery_FDR = if (nrow(discovery_row)) discovery_row$discovery_FDR[[1]] else NA_real_,
      discovery_direction = discovery_direction,
      expansion_logFC = expansion_logfc,
      expansion_PValue = expansion_p,
      expansion_FDR = expansion_fdr,
      expansion_direction = expansion_direction,
      direction_concordant = direction_concordant,
      expansion_nominal_P_lt_0_05 = is.finite(expansion_p) && expansion_p < 0.05,
      expansion_FDR_lt_0_05 = is.finite(expansion_fdr) && expansion_fdr < 0.05,
      source_file = if (!is.null(hit) && nrow(hit)) result_path else NA_character_,
      stringsAsFactors = FALSE
    )
  }
}

coverage <- do.call(rbind, rows)
coverage <- coverage[order(coverage$subtype, match(coverage$gene, panel)), , drop = FALSE]
rownames(coverage) <- NULL
write.csv(
  coverage,
  file.path(out_dir, "23例候选节点覆盖率审计表.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

summary <- do.call(rbind, lapply(split(coverage, coverage$gene), function(x) {
  data.frame(
    gene = x$gene[[1]],
    subtype_records = nrow(x),
    subtypes_with_feature = sum(x$feature_present_in_count_matrix, na.rm = TRUE),
    subtypes_with_reported_effect = sum(is.finite(x$expansion_logFC)),
    concordant_subtypes = sum(x$direction_concordant %in% TRUE, na.rm = TRUE),
    discordant_subtypes = sum(x$direction_concordant %in% FALSE, na.rm = TRUE),
    nominal_positive_subtypes = sum(x$expansion_nominal_P_lt_0_05 %in% TRUE, na.rm = TRUE),
    fdr_positive_subtypes = sum(x$expansion_FDR_lt_0_05 %in% TRUE, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
summary <- summary[match(panel, summary$gene), , drop = FALSE]
write.csv(summary, file.path(out_dir, "23例候选节点覆盖率汇总.csv"), row.names = FALSE, fileEncoding = "UTF-8")

report <- c(
  "# 23例候选节点覆盖率审计",
  "",
  paste0("生成时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "本审计只读取文章1补充表S5和S6中的候选节点，以及23例独立扩大队列已生成的患儿级髓系伪批量结果。它不重新挑选基因，也不把细胞层结果当作患儿级重复。",
  "",
  paste0("候选节点数：", length(panel)),
  paste0("23例髓系状态：", paste(subtypes, collapse = ", ")), 
  "每个状态的重症效应、P值和FDR直接读取既有edgeR准似然输出；覆盖率表示节点是否存在于原始计数矩阵，以及在多少个满足细胞数阈值的患儿级伪批量中有正计数。",
  "",
  "HIF1A、THBS1和CD36的发现方向来自10例发现队列S5；23例结果只用于预设方向复核。方向一致不等于统计学显著，FDR仍按原分析结果保留。",
  "",
  "详细记录：23例候选节点覆盖率审计表.csv；按基因汇总：23例候选节点覆盖率汇总.csv"
)
writeLines(report, file.path(out_dir, "候选节点覆盖率审计说明.md"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "candidate_coverage_sessionInfo.txt"), useBytes = TRUE)
message("Candidate coverage audit completed: ", nrow(coverage), " rows")
