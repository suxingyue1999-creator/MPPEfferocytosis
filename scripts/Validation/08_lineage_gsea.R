#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
})

set.seed(20260721)
workspace <- Sys.getenv("MPP_WORKSPACE", "inputs/workspace")
work_root <- file.path(workspace, "工作文件", "23例单细胞重新分析")
input_root <- file.path(work_root, "results", "07_lineage_pseudobulk")
out_root <- file.path(work_root, "results", "08_lineage_gsea")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

wp_gmt <- system.file(
  "extdata",
  "wikipathways-20180810-gmt-Homo_sapiens.gmt",
  package = "clusterProfiler"
)
if (!file.exists(wp_gmt)) stop("Human WikiPathways GMT is missing: ", wp_gmt)
wp_term2gene <- read.gmt(wp_gmt)
colnames(wp_term2gene) <- c("term", "gene")

all_files <- list.files(
  input_root,
  pattern = "_all_genes\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)
if (!length(all_files)) stop("No lineage pseudobulk DE files found")

safe_gsea_go <- function(gene_list) {
  tryCatch(
    gseGO(
      geneList = gene_list,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      minGSSize = 10,
      maxGSSize = 500,
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      verbose = FALSE,
      seed = TRUE
    ),
    error = function(e) e
  )
}

safe_gsea_wiki <- function(gene_list) {
  tryCatch(
    GSEA(
      geneList = gene_list,
      TERM2GENE = wp_term2gene,
      minGSSize = 10,
      maxGSSize = 500,
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      verbose = FALSE,
      seed = TRUE
    ),
    error = function(e) e
  )
}

empty_result <- function() {
  data.frame()
}

write_enrichment <- function(enrichment, path) {
  result <- if (inherits(enrichment, "error") || is.null(enrichment)) empty_result() else as.data.frame(enrichment)
  write.csv(result, path, row.names = FALSE, fileEncoding = "UTF-8")
  result
}

summary_rows <- list()
for (input_file in all_files) {
  de <- read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("gene", "logFC", "F")
  if (!all(required %in% names(de))) next
  de <- de[is.finite(de$logFC) & is.finite(de$F) & de$F >= 0 & nzchar(de$gene), , drop = FALSE]
  if (nrow(de) < 50L) next

  # Signed sqrt(F) approximates a signed test statistic and preserves both tails.
  de$rank_stat <- sign(de$logFC) * sqrt(de$F)
  de <- de[order(abs(de$rank_stat), decreasing = TRUE), , drop = FALSE]
  de <- de[!duplicated(de$gene), , drop = FALSE]
  mapping <- tryCatch(
    suppressWarnings(bitr(unique(de$gene), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)),
    error = function(e) data.frame(SYMBOL = character(), ENTREZID = character())
  )
  if (!nrow(mapping)) next
  rank_data <- merge(de[, c("gene", "rank_stat"), drop = FALSE], mapping, by.x = "gene", by.y = "SYMBOL", all = FALSE)
  rank_data <- rank_data[order(abs(rank_data$rank_stat), decreasing = TRUE), , drop = FALSE]
  rank_data <- rank_data[!duplicated(rank_data$ENTREZID), , drop = FALSE]
  gene_list <- rank_data$rank_stat
  names(gene_list) <- rank_data$ENTREZID
  gene_list <- sort(gene_list, decreasing = TRUE)
  if (length(gene_list) < 20L) next

  rel <- sub(paste0("^", normalizePath(input_root, mustWork = TRUE), "/?"), "", normalizePath(input_file, mustWork = TRUE))
  parts <- strsplit(rel, .Platform$file.sep, fixed = TRUE)[[1]]
  lineage <- parts[1]
  subtype <- parts[2]
  contrast <- sub("_all_genes\\.csv$", "", parts[length(parts)])
  result_dir <- file.path(out_root, lineage, subtype, contrast)
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

  go <- safe_gsea_go(gene_list)
  wiki <- safe_gsea_wiki(gene_list)
  go_df <- write_enrichment(go, file.path(result_dir, "GO_BP_GSEA.csv"))
  wiki_df <- write_enrichment(wiki, file.path(result_dir, "WikiPathways_GSEA.csv"))

  significant <- de[de$FDR < 0.05 & abs(de$logFC) >= 1, "gene", drop = TRUE]
  sig_mapping <- mapping[mapping$SYMBOL %in% significant, , drop = FALSE]
  ora <- if (length(unique(sig_mapping$ENTREZID)) >= 5L) {
    tryCatch(
      enrichGO(
        gene = unique(sig_mapping$ENTREZID),
        universe = unique(mapping$ENTREZID),
        OrgDb = org.Hs.eg.db,
        keyType = "ENTREZID",
        ont = "BP",
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        readable = TRUE
      ),
      error = function(e) NULL
    )
  } else NULL
  ora_df <- write_enrichment(ora, file.path(result_dir, "GO_BP_ORA_significant.csv"))

  if (nrow(go_df)) {
    pdf(file.path(result_dir, "GO_BP_GSEA_dotplot.pdf"), width = 10, height = 7)
    print(dotplot(go, showCategory = min(15, nrow(go_df))) + ggtitle(paste(lineage, subtype, contrast, "GO BP")))
    dev.off()
  }
  if (nrow(wiki_df)) {
    pdf(file.path(result_dir, "WikiPathways_GSEA_dotplot.pdf"), width = 10, height = 7)
    print(dotplot(wiki, showCategory = min(15, nrow(wiki_df))) + ggtitle(paste(lineage, subtype, contrast, "WikiPathways")))
    dev.off()
  }

  summary_rows[[length(summary_rows) + 1L]] <- data.frame(
    lineage = lineage,
    subtype = subtype,
    contrast = contrast,
    input_genes = nrow(de),
    mapped_entrez = length(gene_list),
    significant_genes = length(significant),
    GO_BP_terms = nrow(go_df),
    GO_BP_FDR_lt_0_05 = if (nrow(go_df)) sum(go_df$p.adjust < 0.05, na.rm = TRUE) else 0L,
    WikiPathways_terms = nrow(wiki_df),
    WikiPathways_FDR_lt_0_05 = if (nrow(wiki_df)) sum(wiki_df$p.adjust < 0.05, na.rm = TRUE) else 0L,
    GO_BP_ORA_terms = nrow(ora_df),
    status = "completed",
    stringsAsFactors = FALSE
  )
}

summary <- do.call(rbind, summary_rows)
write.csv(summary, file.path(out_root, "lineage_gsea_summary.csv"), row.names = FALSE, fileEncoding = "UTF-8")

report_lines <- c(
  "# 四大谱系通路分析",
  "",
  paste0("- 运行时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("- 输入：", nrow(summary), " 个完成的亚群×临床对比 edgeR 全基因结果。"),
  "- GSEA 排名：`sign(logFC) × sqrt(F)`，使用所有通过 edgeR 表达过滤的基因，不只使用显著基因。",
  "- GO：使用本机 `clusterProfiler` + `org.Hs.eg.db` 做离线 GO Biological Process GSEA；另对达到 FDR<0.05 且 |log2FC|>=1 的基因做 GO-BP ORA。",
  paste0("- WikiPathways：使用本机 GMT `", wp_gmt, "`；该资源版本为 2018-08-10，解释时应标注资源版本。"),
  "- 未使用 msigdbr 或 ReactomePA；没有上传表达数据，所有排名和富集计算在本机完成。",
  "- 中性粒细胞模块标签全部为 inflammatory，因此其通路结果只能解释为该标签下的 RNA 状态，不等同于多个可靠中性粒细胞亚型。",
  "",
  "`lineage_gsea_summary.csv` 汇总每个比较的输入基因、Entrez 映射率、GO/WikiPathways 富集条目及 FDR<0.05 数量；每个比较的完整结果在对应子目录。"
)
writeLines(report_lines, file.path(out_root, "四大谱系通路分析说明.md"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(out_root, "sessionInfo.txt"))
message("Lineage GO and WikiPathways GSEA completed")
