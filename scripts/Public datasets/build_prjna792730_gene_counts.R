#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("Usage: build_prjna792730_gene_counts.R <gencode.gtf.gz> <quant_dir> <count_output.tsv> <qc_output.tsv>")
}

gtf_path <- normalizePath(args[[1L]], mustWork = TRUE)
quant_dir <- normalizePath(args[[2L]], mustWork = TRUE)
count_output <- args[[3L]]
qc_output <- args[[4L]]

sample_info <- data.table(
  run = c("SRR17345233", "SRR17345232", "SRR17345231",
          "SRR17345230", "SRR17345229", "SRR17345228",
          "SRR17345227", "SRR17345226", "SRR17345225"),
  sample = c("C1-0h", "C2-0h", "C3-0h", "T1-12h", "T2-12h", "T3-12h",
             "T1-24h", "T2-24h", "T3-24h"),
  timepoint = c("0h", "0h", "0h", "12h", "12h", "12h", "24h", "24h", "24h")
)

gtf <- fread(
  gtf_path,
  sep = "\t",
  header = FALSE,
  quote = "",
  select = c(3L, 9L),
  col.names = c("feature", "attributes")
)
gtf <- gtf[feature == "transcript"]
extract_attribute <- function(x, key) {
  pattern <- paste0("(?:^|;[[:space:]]*)", key, " \\\"([^\\\"]+)\\\"")
  sub(paste0(".*", pattern, ".*"), "\\1", x, perl = TRUE)
}
gtf[, transcript_id := extract_attribute(attributes, "transcript_id")]
gtf[, gene_id := extract_attribute(attributes, "gene_id")]
gtf[, gene_name := extract_attribute(attributes, "gene_name")]
# In GENCODE v44 the transcript_id field already carries the Ensembl version
# suffix, matching the transcript identifier in the FASTA header used by Salmon.
gtf[, transcript := transcript_id]
tx2gene <- unique(gtf[, .(transcript, gene_id, gene_name)])
if (anyNA(tx2gene$transcript) || anyNA(tx2gene$gene_id) || anyNA(tx2gene$gene_name) ||
    any(!nzchar(tx2gene$transcript)) || any(!nzchar(tx2gene$gene_id)) || any(!nzchar(tx2gene$gene_name))) {
  stop("GENCODE transcript-to-gene mapping contains missing or empty identifiers.")
}
if (anyDuplicated(tx2gene$transcript)) {
  stop("GENCODE transcript identifiers are not uniquely mapped to genes.")
}

gene_counts <- NULL
qc_rows <- vector("list", nrow(sample_info))
for (i in seq_len(nrow(sample_info))) {
  sample <- sample_info$sample[[i]]
  quant_path <- file.path(quant_dir, sample, "quant.sf")
  meta_path <- file.path(quant_dir, sample, "aux_info", "meta_info.json")
  if (!file.exists(quant_path) || !file.exists(meta_path)) {
    stop("Missing completed Salmon output for ", sample)
  }
  quant <- fread(quant_path, select = c("Name", "TPM", "NumReads"))
  quant[, transcript := sub("\\|.*$", "", Name)]
  quant <- merge(quant, tx2gene, by = "transcript", all.x = TRUE, sort = FALSE)
  mapping_rate <- mean(!is.na(quant$gene_name))
  if (mapping_rate < 0.99) {
    stop(sprintf("Only %.2f%% of Salmon targets mapped to GENCODE genes for %s.",
                 100 * mapping_rate, sample))
  }
  by_gene <- quant[!is.na(gene_name), .(
    estimated_count = sum(NumReads),
    TPM = sum(TPM)
  ), by = .(gene_id, gene_name)]
  setnames(by_gene, "estimated_count", sample)
  if (is.null(gene_counts)) {
    gene_counts <- by_gene[, c("gene_id", "gene_name", sample), with = FALSE]
  } else {
    gene_counts <- merge(
      gene_counts,
      by_gene[, c("gene_id", "gene_name", sample), with = FALSE],
      by = c("gene_id", "gene_name"),
      all = TRUE,
      sort = FALSE
    )
  }

  meta <- fromJSON(meta_path)
  qc_rows[[i]] <- data.table(
    run = sample_info$run[[i]],
    sample = sample,
    timepoint = sample_info$timepoint[[i]],
    num_processed = as.numeric(meta$num_processed),
    num_mapped = as.numeric(meta$num_mapped),
    mapping_rate = as.numeric(meta$percent_mapped) / 100,
    inferred_library_type = as.character(meta$library_types),
    salmon_version = as.character(meta$salmon_version),
    gencode_target_mapping_rate = mapping_rate
  )
}

for (column in sample_info$sample) {
  set(gene_counts, which(is.na(gene_counts[[column]])), column, 0)
}
setorder(gene_counts, gene_name, gene_id)
fwrite(gene_counts, count_output, sep = "\t", quote = FALSE, na = "NA")
fwrite(rbindlist(qc_rows), qc_output, sep = "\t", quote = FALSE, na = "NA")
