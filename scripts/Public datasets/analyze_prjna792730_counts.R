#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop("Usage: analyze_prjna792730_counts.R <count_matrix.tsv[.gz]> <output_dir>")
}

count_path <- normalizePath(args[[1L]], mustWork = TRUE)
out_dir <- args[[2L]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

candidate_genes <- c(
  "HIF1A", "THBS1", "CD36", "TGFB1", "MAPK1", "MAPK3", "SMAD2", "SMAD4",
  "MERTK", "AXL", "C1QA", "C1QB", "C1QC", "APOE", "TREM2", "CTSD",
  "LAMP1", "NPC2", "LIPA"
)

sample_info <- data.table(
  sample = c("C1-0h", "C2-0h", "C3-0h", "T1-12h", "T2-12h", "T3-12h",
             "T1-24h", "T2-24h", "T3-24h"),
  timepoint = factor(c("0h", "0h", "0h", "12h", "12h", "12h", "24h", "24h", "24h"),
                     levels = c("0h", "12h", "24h")),
  replicate = c("C1", "C2", "C3", "T1", "T2", "T3", "T1", "T2", "T3")
)

counts_dt <- fread(count_path)
if (ncol(counts_dt) < 10L) {
  stop("The count matrix must contain one gene identifier column and all nine PRJNA792730 samples.")
}

sample_names <- setdiff(names(counts_dt), c("gene_id", "gene_name"))
missing_samples <- setdiff(sample_info$sample, sample_names)
if (length(missing_samples)) {
  stop("Missing expected samples: ", paste(missing_samples, collapse = ", "))
}
extra_samples <- setdiff(sample_names, sample_info$sample)
if (length(extra_samples)) {
  warning("Ignoring columns not present in the prespecified nine-sample design: ",
          paste(extra_samples, collapse = ", "))
}

if (!all(c("gene_id", "gene_name") %in% names(counts_dt))) {
  stop("The count matrix must retain gene_id and gene_name annotation columns.")
}
counts_dt <- counts_dt[, c("gene_id", "gene_name", sample_info$sample), with = FALSE]
gene_ids <- as.character(counts_dt$gene_id)
gene_names <- as.character(counts_dt$gene_name)
if (anyNA(gene_ids) || any(!nzchar(gene_ids)) || anyDuplicated(gene_ids)) {
  stop("gene_id must be non-missing and unique in the count matrix.")
}
counts <- do.call(cbind, lapply(counts_dt[, sample_info$sample, with = FALSE], as.numeric))
if (anyNA(counts) || any(counts < 0)) {
  stop("The input matrix must contain non-negative estimated counts without missing values.")
}
rownames(counts) <- gene_ids
colnames(counts) <- sample_info$sample
storage.mode(counts) <- "double"

gene_name_by_id <- setNames(gene_names, gene_ids)
candidate_map <- data.table(gene = gene_names, gene_id = gene_ids)[gene %in% candidate_genes]
candidate_map <- candidate_map[, .(gene_ids = list(gene_id), n_ids = .N), by = gene]
ambiguous_candidates <- candidate_map[n_ids > 1L, gene]
if (length(ambiguous_candidates)) {
  stop("Candidate gene symbols map to multiple gene_id values: ", paste(ambiguous_candidates, collapse = ", "))
}
candidate_id <- setNames(vapply(candidate_genes, function(g) {
  ids <- candidate_map[gene == g, unlist(gene_ids)]
  if (length(ids) == 1L) ids else NA_character_
}, character(1)), candidate_genes)

## The public study is a bulk THP-1 macrophage-like time series. The design is
## prespecified as timepoint only; no pediatric severity label is available.
design <- model.matrix(~ 0 + timepoint, data = sample_info)
colnames(design) <- c("t0h", "t12h", "t24h")
y_all <- DGEList(counts = counts, samples = sample_info)
keep <- filterByExpr(y_all, design = design)
if (!any(keep)) stop("No genes passed edgeR filterByExpr for the nine-sample matrix.")
y <- calcNormFactors(y_all[keep, , keep.lib.sizes = FALSE])
fit <- glmQLFit(y, design, robust = TRUE)

contrasts <- list(
  `12h_vs_0h` = makeContrasts(t12h - t0h, levels = design),
  `24h_vs_0h` = makeContrasts(t24h - t0h, levels = design),
  `24h_vs_12h` = makeContrasts(t24h - t12h, levels = design)
)

all_results <- rbindlist(lapply(names(contrasts), function(label) {
  qlf <- glmQLFTest(fit, contrast = contrasts[[label]])
  tab <- as.data.table(topTags(qlf, n = Inf, sort.by = "none")$table, keep.rownames = "gene_id")
  tab[, gene_name := unname(gene_name_by_id[gene_id])]
  tab[, `:=`(dataset = "PRJNA792730", comparison = label,
             comparison_role = fifelse(label == "24h_vs_0h", "primary_publication_aligned", "secondary_time_response"),
             tested_genes = nrow(tab))]
  setcolorder(tab, c("dataset", "comparison", "comparison_role", "gene_id", "gene_name", "logFC", "logCPM", "F", "PValue", "FDR", "tested_genes"))
  tab
}), use.names = TRUE, fill = TRUE)
fwrite(all_results, file.path(out_dir, "PRJNA792730_all_gene_results.csv"), na = "NA")

primary_result <- all_results[comparison == "24h_vs_0h"]
publication_concordance <- data.table(
  dataset = "PRJNA792730",
  comparison = "24h_vs_0h",
  threshold = "absolute log2 fold change > 1 and whole-transcriptome FDR < 0.05",
  reanalysis_up = primary_result[logFC > 1 & FDR < 0.05, .N],
  reanalysis_down = primary_result[logFC < -1 & FDR < 0.05, .N],
  reanalysis_total = primary_result[abs(logFC) > 1 & FDR < 0.05, .N],
  publication_reported_up = 927L,
  publication_reported_down = 1162L,
  publication_reported_total = 2089L,
  publication_source = "Ma et al., International Journal of Molecular Sciences, 2023, DOI 10.3390/ijms24108588"
)
fwrite(
  publication_concordance,
  file.path(out_dir, "PRJNA792730_publication_concordance_audit.csv"),
  na = "NA"
)

logcpm <- cpm(y, log = TRUE, prior.count = 1)
candidate_results <- rbindlist(lapply(names(contrasts), function(label) {
  qlf <- glmQLFTest(fit, contrast = contrasts[[label]])
  tab <- topTags(qlf, n = Inf, sort.by = "none")$table
  out <- data.table(
    dataset = "PRJNA792730",
    comparison = label,
    comparison_role = if (label == "24h_vs_0h") "primary_publication_aligned" else "secondary_time_response",
    gene = candidate_genes,
    detected_n = vapply(candidate_genes, function(g) {
      id <- candidate_id[[g]]
      if (!is.na(id) && id %in% rownames(counts)) sum(counts[id, ] > 0) else 0L
    }, integer(1)),
    detected_rate = vapply(candidate_genes, function(g) {
      id <- candidate_id[[g]]
      if (!is.na(id) && id %in% rownames(counts)) mean(counts[id, ] > 0) else 0
    }, numeric(1)),
    mean_logCPM_0h = vapply(candidate_genes, function(g) {
      id <- candidate_id[[g]]
      if (!is.na(id) && id %in% rownames(logcpm)) mean(logcpm[id, sample_info$timepoint == "0h"]) else NA_real_
    }, numeric(1)),
    mean_logCPM_12h = vapply(candidate_genes, function(g) {
      id <- candidate_id[[g]]
      if (!is.na(id) && id %in% rownames(logcpm)) mean(logcpm[id, sample_info$timepoint == "12h"]) else NA_real_
    }, numeric(1)),
    mean_logCPM_24h = vapply(candidate_genes, function(g) {
      id <- candidate_id[[g]]
      if (!is.na(id) && id %in% rownames(logcpm)) mean(logcpm[id, sample_info$timepoint == "24h"]) else NA_real_
    }, numeric(1)),
    logFC = vapply(candidate_genes, function(g) {
      id <- candidate_id[[g]]
      if (!is.na(id) && id %in% rownames(tab)) unname(tab[id, "logFC"]) else NA_real_
    }, numeric(1)),
    P_value = vapply(candidate_genes, function(g) {
      id <- candidate_id[[g]]
      if (!is.na(id) && id %in% rownames(tab)) unname(tab[id, "PValue"]) else NA_real_
    }, numeric(1)),
    FDR = vapply(candidate_genes, function(g) {
      id <- candidate_id[[g]]
      if (!is.na(id) && id %in% rownames(tab)) unname(tab[id, "FDR"]) else NA_real_
    }, numeric(1))
  )
  out[, candidate_family_tested_n := sum(!is.na(P_value))]
  out[, FDR_candidate_tested := p.adjust(P_value, method = "BH")]
  out
}), use.names = TRUE, fill = TRUE)
fwrite(candidate_results, file.path(out_dir, "PRJNA792730_candidate_nodes.csv"), na = "NA")
fwrite(sample_info, file.path(out_dir, "PRJNA792730_sample_metadata.csv"), na = "NA")

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
cat("Completed prespecified PRJNA792730 time-course analysis.\n")
