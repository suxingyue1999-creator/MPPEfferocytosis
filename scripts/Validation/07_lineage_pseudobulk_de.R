#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(edgeR)
  library(ggplot2)
})

set.seed(20260721)
options(future.globals.maxSize = 40 * 1024^3)
future::plan("sequential")

workspace <- Sys.getenv("MPP_WORKSPACE", "inputs/workspace")
work_root <- file.path(workspace, "工作文件", "23例单细胞重新分析")
input_root <- file.path(work_root, "results", "02_rna_reanalysis", "subclusters")
out_root <- file.path(work_root, "results", "07_lineage_pseudobulk")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

lineage_stems <- c("Macrophage_Monocyte", "T_cells", "B_cell", "Neutrophils")
required_meta <- c(
  "sample_id", "severity", "allergy_group", "plastic_group", "batch_wave",
  "subcluster_id", "subtype_label"
)
clinical_columns <- c(
  "sample_id", "severity", "allergy_group", "plastic_group", "batch_wave"
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

read_lineage <- function(stem) {
  path <- file.path(input_root, stem, paste0(stem, "_subclustered.rds"))
  if (!file.exists(path)) stop("Missing lineage object: ", path)
  object <- readRDS(path)
  if (!inherits(object, "Seurat")) stop("Not a Seurat object: ", path)
  DefaultAssay(object) <- "RNA"
  missing <- setdiff(required_meta, colnames(object[[]]))
  if (length(missing)) stop(stem, " missing metadata: ", paste(missing, collapse = ", "))
  object <- DietSeurat(
    object,
    assays = "RNA",
    layers = "counts",
    dimreducs = NULL,
    graphs = NULL
  )
  object
}

prepare_metadata <- function(meta) {
  meta <- meta[, unique(c(required_meta, clinical_columns)), drop = FALSE]
  meta$sample_id <- as.character(meta$sample_id)
  meta$severity <- factor(as.character(meta$severity), levels = c("轻症", "重症"))
  meta$allergy_group <- factor(as.character(meta$allergy_group), levels = c("否", "是"))
  meta$plastic_group <- factor(as.character(meta$plastic_group), levels = c("否", "是"))
  meta$batch_wave <- factor(as.character(meta$batch_wave), levels = c("2023_Dec", "2024_Spring"))
  meta$subtype_label <- as.character(meta$subtype_label)
  meta
}

make_design <- function(meta) {
  model.matrix(~ batch_wave + severity + allergy_group + plastic_group, data = meta)
}

write_skip_rows <- function(celltype, contrast_names, samples, min_cells, rank, columns, status) {
  do.call(rbind, lapply(contrast_names, function(contrast_name) {
    data.frame(
      subtype = celltype,
      contrast = contrast_labels[[contrast_name]],
      samples_used = samples,
      min_cells_per_pseudobulk = min_cells,
      design_rank = rank,
      design_columns = columns,
      genes_tested = NA_integer_,
      FDR_lt_0_05 = NA_integer_,
      FDR_lt_0_05_abs_logFC_ge_1 = NA_integer_,
      status = status,
      stringsAsFactors = FALSE
    )
  }))
}

all_de_summary <- list()
all_composition_summary <- list()
summary_index <- 0L
composition_index <- 0L

for (stem in lineage_stems) {
  message("Lineage pseudobulk: ", stem)
  lineage_out <- file.path(out_root, stem)
  dir.create(lineage_out, recursive = TRUE, showWarnings = FALSE)
  object <- read_lineage(stem)
  meta <- prepare_metadata(object[[]])
  counts <- LayerData(object[["RNA"]], layer = "counts")
  if (!identical(colnames(counts), rownames(meta))) {
    meta <- meta[colnames(counts), , drop = FALSE]
  }
  if (!identical(colnames(counts), rownames(meta))) {
    stop(stem, ": counts and metadata cell order do not match")
  }

  sample_info <- unique(meta[, clinical_columns, drop = FALSE])
  if (anyDuplicated(sample_info$sample_id)) stop(stem, ": sample metadata is not unique")
  sample_info <- sample_info[order(sample_info$sample_id), , drop = FALSE]
  subtypes <- sort(unique(meta$subtype_label))
  group_id <- paste(meta$sample_id, meta$subtype_label, sep = "||")
  group_factor <- factor(group_id, levels = unique(group_id))
  design_cells <- sparse.model.matrix(~ 0 + group_factor)
  colnames(design_cells) <- levels(group_factor)
  pseudobulk_counts <- counts %*% design_cells
  rm(counts, design_cells)
  invisible(gc(verbose = FALSE))

  group_parts <- strsplit(colnames(pseudobulk_counts), "\\|\\|")
  pseudobulk_meta <- data.frame(
    pseudobulk_id = colnames(pseudobulk_counts),
    sample_id = vapply(group_parts, `[[`, character(1), 1L),
    subtype = vapply(group_parts, `[[`, character(1), 2L),
    stringsAsFactors = FALSE
  )
  group_sizes <- table(group_factor)
  pseudobulk_meta$cell_count <- as.integer(group_sizes[pseudobulk_meta$pseudobulk_id])
  pseudobulk_meta <- cbind(
    pseudobulk_meta,
    sample_info[match(pseudobulk_meta$sample_id, sample_info$sample_id), clinical_columns[-1], drop = FALSE]
  )
  write.csv(
    pseudobulk_meta,
    file.path(lineage_out, "pseudobulk_sample_metadata.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
  saveRDS(
    pseudobulk_counts,
    file.path(lineage_out, "pseudobulk_raw_counts.rds"),
    compress = "gzip"
  )

  lineage_counts <- as.data.frame.matrix(table(meta$sample_id, meta$subtype_label))
  lineage_counts <- lineage_counts[match(sample_info$sample_id, rownames(lineage_counts)), , drop = FALSE]
  lineage_counts$sample_id <- rownames(lineage_counts)
  rownames(lineage_counts) <- NULL
  lineage_counts <- lineage_counts[, c("sample_id", setdiff(names(lineage_counts), "sample_id")), drop = FALSE]
  lineage_count_matrix <- as.matrix(lineage_counts[, setdiff(names(lineage_counts), "sample_id"), drop = FALSE])
  lineage_total <- rowSums(lineage_count_matrix)
  composition <- do.call(rbind, lapply(seq_len(nrow(lineage_counts)), function(i) {
    data.frame(
      sample_id = lineage_counts$sample_id[i],
      subtype = colnames(lineage_count_matrix),
      cell_count = as.integer(lineage_count_matrix[i, ]),
      lineage_total = as.integer(lineage_total[i]),
      stringsAsFactors = FALSE
    )
  }))
  composition$fraction <- composition$cell_count / composition$lineage_total
  composition <- merge(composition, sample_info, by = "sample_id", all.x = TRUE, sort = FALSE)
  composition <- composition[order(match(composition$sample_id, sample_info$sample_id), composition$subtype), , drop = FALSE]
  write.csv(composition, file.path(lineage_out, "subtype_composition_by_sample.csv"), row.names = FALSE, fileEncoding = "UTF-8")

  for (subtype in subtypes) {
    keep_columns <- which(pseudobulk_meta$subtype == subtype & pseudobulk_meta$cell_count >= 20L)
    current_meta <- pseudobulk_meta[keep_columns, , drop = FALSE]
    current_counts <- pseudobulk_counts[, keep_columns, drop = FALSE]
    if (!nrow(current_meta)) {
      summary_index <- summary_index + 1L
      all_de_summary[[summary_index]] <- write_skip_rows(subtype, names(contrast_map), 0L, NA_integer_, NA_integer_, NA_integer_, "skipped_no_pseudobulk")
      next
    }
    current_meta$batch_wave <- factor(as.character(current_meta$batch_wave), levels = c("2023_Dec", "2024_Spring"))
    current_meta$severity <- factor(as.character(current_meta$severity), levels = c("轻症", "重症"))
    current_meta$allergy_group <- factor(as.character(current_meta$allergy_group), levels = c("否", "是"))
    current_meta$plastic_group <- factor(as.character(current_meta$plastic_group), levels = c("否", "是"))
    design <- make_design(current_meta)
    design_rank <- qr(design)$rank
    design_columns <- ncol(design)
    min_cells <- min(current_meta$cell_count)
    subtype_dir <- file.path(lineage_out, gsub("[^A-Za-z0-9_]+", "_", subtype))
    dir.create(subtype_dir, recursive = TRUE, showWarnings = FALSE)
    write.csv(
      data.frame(sample_id = current_meta$sample_id, design, check.names = FALSE),
      file.path(subtype_dir, "design_matrix.csv"), row.names = FALSE
    )
    if (nrow(current_meta) < 8L || design_rank != design_columns) {
      summary_index <- summary_index + 1L
      all_de_summary[[summary_index]] <- write_skip_rows(subtype, names(contrast_map), nrow(current_meta), min_cells, design_rank, design_columns, "skipped_insufficient_samples_or_rank")
      next
    }

    dge <- DGEList(counts = current_counts)
    keep_genes <- filterByExpr(dge, design = design)
    dge <- dge[keep_genes, , keep.lib.sizes = FALSE]
    if (nrow(dge) < 10L) {
      summary_index <- summary_index + 1L
      all_de_summary[[summary_index]] <- write_skip_rows(subtype, names(contrast_map), nrow(current_meta), min_cells, design_rank, design_columns, "skipped_too_few_genes_after_filter")
      next
    }
    dge <- calcNormFactors(dge, method = "TMM")
    dge <- estimateDisp(dge, design = design, robust = TRUE)
    fit <- glmQLFit(dge, design = design, robust = TRUE)

    for (contrast_name in names(contrast_map)) {
      coefficient <- contrast_map[[contrast_name]]
      label <- contrast_labels[[contrast_name]]
      if (!coefficient %in% colnames(design)) {
        summary_index <- summary_index + 1L
        all_de_summary[[summary_index]] <- write_skip_rows(subtype, contrast_name, nrow(current_meta), min_cells, design_rank, design_columns, "skipped_coefficient_absent")
        next
      }
      test <- glmQLFTest(fit, coef = coefficient)
      result <- topTags(test, n = Inf, sort.by = "PValue")$table
      result$gene <- rownames(result)
      rownames(result) <- NULL
      result <- result[, c("gene", setdiff(names(result), "gene")), drop = FALSE]
      result$significant_FDR_0_05 <- result$FDR < 0.05
      result$significant_FDR_0_05_abs_logFC_1 <- result$FDR < 0.05 & abs(result$logFC) >= 1
      stem_name <- paste0(contrast_name, "_", gsub("[^A-Za-z0-9_]+", "_", label))
      write.csv(result, file.path(subtype_dir, paste0(stem_name, "_all_genes.csv")), row.names = FALSE, fileEncoding = "UTF-8")
      write.csv(result[result$significant_FDR_0_05_abs_logFC_1, , drop = FALSE], file.path(subtype_dir, paste0(stem_name, "_FDR005_logFC1.csv")), row.names = FALSE, fileEncoding = "UTF-8")
      summary_index <- summary_index + 1L
      all_de_summary[[summary_index]] <- data.frame(
        subtype = subtype,
        contrast = label,
        samples_used = nrow(current_meta),
        min_cells_per_pseudobulk = min_cells,
        design_rank = design_rank,
        design_columns = design_columns,
        genes_tested = nrow(result),
        FDR_lt_0_05 = sum(result$FDR < 0.05, na.rm = TRUE),
        FDR_lt_0_05_abs_logFC_ge_1 = sum(result$significant_FDR_0_05_abs_logFC_1, na.rm = TRUE),
        status = "completed",
        stringsAsFactors = FALSE
      )
    }
  }

  composition_results <- list()
  for (subtype in subtypes) {
    current <- composition[composition$subtype == subtype, , drop = FALSE]
    if (nrow(current) < 8L || any(current$lineage_total <= 0)) next
    design <- make_design(current)
    rank <- qr(design)$rank
    if (rank != ncol(design)) next
    fit <- tryCatch(
      glm(cbind(cell_count, lineage_total - cell_count) ~ batch_wave + severity + allergy_group + plastic_group, family = quasibinomial(), data = current),
      error = function(e) NULL
    )
    if (is.null(fit)) next
    coefficient_table <- summary(fit)$coefficients
    for (contrast_name in names(contrast_map)) {
      coefficient <- contrast_map[[contrast_name]]
      if (!coefficient %in% rownames(coefficient_table)) next
      composition_index <- composition_index + 1L
      composition_results[[composition_index]] <- data.frame(
        lineage = stem,
        subtype = subtype,
        contrast = contrast_labels[[contrast_name]],
        samples_used = nrow(current),
        estimate_log_odds = coefficient_table[coefficient, "Estimate"],
        odds_ratio = exp(coefficient_table[coefficient, "Estimate"]),
        p_value = coefficient_table[coefficient, "Pr(>|t|)"],
        status = "exploratory_quasibinomial",
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(composition_results)) {
    composition_summary <- do.call(rbind, composition_results)
    composition_summary$FDR <- p.adjust(composition_summary$p_value, method = "BH")
    write.csv(composition_summary, file.path(lineage_out, "subtype_composition_quasibinomial.csv"), row.names = FALSE, fileEncoding = "UTF-8")
    all_composition_summary[[stem]] <- composition_summary
  }
  for (variable in c("severity", "allergy_group", "plastic_group")) {
    plot <- ggplot(composition, aes(x = .data[[variable]], y = fraction, color = .data[[variable]])) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.12, height = 0, alpha = 0.75, size = 1.2) +
      facet_wrap(~ subtype, scales = "free_y") +
      theme_bw(base_size = 11) +
      labs(x = variable, y = "Sample fraction", title = paste(stem, "subtype composition by", variable)) +
      theme(legend.position = "none")
    ggsave(file.path(lineage_out, paste0("subtype_composition_", variable, ".pdf")), plot, width = 12, height = 7, useDingbats = FALSE)
    ggsave(file.path(lineage_out, paste0("subtype_composition_", variable, ".png")), plot, width = 12, height = 7, dpi = 180)
  }
  rm(object, meta, pseudobulk_counts, pseudobulk_meta, composition)
  invisible(gc(verbose = FALSE))
}

de_summary <- do.call(rbind, all_de_summary)
write.csv(de_summary, file.path(out_root, "lineage_pseudobulk_DE_summary.csv"), row.names = FALSE, fileEncoding = "UTF-8")
if (length(all_composition_summary)) {
  composition_summary <- do.call(rbind, all_composition_summary)
  rownames(composition_summary) <- NULL
  write.csv(composition_summary, file.path(out_root, "lineage_subtype_composition_summary.csv"), row.names = FALSE, fileEncoding = "UTF-8")
}

report_lines <- c(
  "# 四大谱系样本级 pseudobulk 和组成分析",
  "",
  paste0("- 运行时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "- 输入：四个谱系的 RNA 子群 RDS；表达输入为 RNA 原始 counts。",
  "- pseudobulk 单位：样本 × RNA 子群；每个组合至少 20 个细胞才进入 DE 模型。",
  "- DE 模型：`~ batch_wave + severity + allergy_group + plastic_group`；edgeR quasi-likelihood、TMM、robust dispersion。",
  "- 组成模型：样本级 `quasibinomial`，以某亚群细胞数对该谱系其余细胞数建模；这是探索性组成检验，不把单细胞当独立重复。",
  "- 对每个谱系均保存原始 pseudobulk counts、样本元数据、全基因 DE 结果、显著结果、组成表和组间比例图。",
  "- 中性粒细胞的模块标签目前全部为 `inflammatory`；该标签不应被解读为已经得到多个可靠中性粒细胞亚型。",
  "",
  "`lineage_pseudobulk_DE_summary.csv` 是 DE 总表；通路分析脚本将使用每个完成模型的全基因排名进行 GSEA。"
)
writeLines(report_lines, file.path(out_root, "四大谱系pseudobulk分析说明.md"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(out_root, "sessionInfo.txt"))
message("Lineage pseudobulk DE and composition analysis completed")
