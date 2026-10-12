#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

base_dir <- Sys.getenv("MPP_PUBLIC_ROOT", "inputs/public")
out_dir <- file.path(base_dir, "results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

candidate_genes <- c(
  "HIF1A", "THBS1", "CD36", "TGFB1", "MAPK1", "MAPK3", "SMAD2", "SMAD4",
  "MERTK", "AXL", "C1QA", "C1QB", "C1QC", "APOE", "TREM2", "CTSD",
  "LAMP1", "NPC2", "LIPA"
)

write_table <- function(x, path) {
  fwrite(x, path, na = "NA", quote = TRUE)
}

## GSE274269: 174 M. pneumoniae-positive BALF total-RNA samples.
## The GEO annotation contains only P1 genotype; no severity or cell-type labels.
gse274_path <- file.path(base_dir, "GSE274269_allsample_featurecount.txt.gz")
## Read the sample identifiers as text so long identifiers are not coerced to doubles.
gse274_header <- strsplit(readLines(gzfile(gse274_path), n = 1L), "\t", fixed = TRUE)[[1L]]
gse274_counts <- fread(
  cmd = paste("gzip -cd", shQuote(gse274_path)),
  skip = 1L,
  data.table = TRUE,
  header = FALSE
)
setnames(gse274_counts, gse274_header)
stopifnot(ncol(gse274_counts) == 175L, gse274_header[1L] == "Geneid")
gene_ids <- gse274_counts[[1L]]
count_mat <- do.call(cbind, lapply(gse274_counts[, -1L, with = FALSE], as.numeric))
colnames(count_mat) <- names(gse274_counts)[-1L]
rownames(count_mat) <- gene_ids
storage.mode(count_mat) <- "numeric"

gse274_series_lines <- readLines(gzfile(file.path(base_dir, "GSE274269_series_matrix.txt.gz")))
gse274_genotype_line <- gse274_series_lines[
  grepl("^!Sample_characteristics_ch1\\t.*genotype:", gse274_series_lines)
]
stopifnot(length(gse274_genotype_line) == 1L)
gse274_genotype <- sub(
  "^genotype: ", "",
  gsub('^"|"$', "", strsplit(gse274_genotype_line, "\\t", fixed = FALSE)[[1L]][-1L])
)
gse274_meta <- data.table(sample_column = colnames(count_mat), genotype = gse274_genotype)
stopifnot(nrow(gse274_meta) == ncol(count_mat),
          sum(gse274_meta$genotype == "P1-1") == 108L,
          sum(gse274_meta$genotype == "P1-2") == 66L)

## Sum duplicated symbols before the candidate-gene lookup.
unique_genes <- !duplicated(rownames(count_mat))
if (!all(unique_genes)) {
  count_dt <- as.data.table(count_mat, keep.rownames = "gene")
  count_dt <- count_dt[, lapply(.SD, sum), by = gene]
  count_mat <- as.matrix(count_dt[, -1L, with = FALSE])
  rownames(count_mat) <- count_dt$gene
  storage.mode(count_mat) <- "numeric"
}

design274 <- model.matrix(~ genotype, data = gse274_meta)
y274_all <- DGEList(counts = count_mat, samples = gse274_meta)
y274_all <- calcNormFactors(y274_all)
logcpm274 <- cpm(y274_all, log = TRUE, prior.count = 1)
keep274 <- filterByExpr(y274_all, design = design274)
y274 <- y274_all[keep274, , keep.lib.sizes = FALSE]
y274 <- calcNormFactors(y274)
fit274 <- glmQLFit(y274, design274, robust = TRUE)
qlf274 <- glmQLFTest(fit274, coef = "genotypeP1-2")
tab274 <- topTags(qlf274, n = Inf, sort.by = "none")$table
tab274$gene <- rownames(tab274)

idx274 <- match(candidate_genes, rownames(logcpm274))
res274 <- data.table(
  dataset = "GSE274269",
  comparison = "P1-2_vs_P1-1_exploratory",
  gene = candidate_genes,
  detected_n = vapply(seq_along(candidate_genes), function(i) {
    if (is.na(idx274[i])) return(0L)
    sum(count_mat[idx274[i], ] > 0)
  }, integer(1)),
  detected_rate = vapply(seq_along(candidate_genes), function(i) {
    if (is.na(idx274[i])) return(0)
    mean(count_mat[idx274[i], ] > 0)
  }, numeric(1)),
  mean_logCPM_P1_1 = vapply(seq_along(candidate_genes), function(i) {
    if (is.na(idx274[i])) return(NA_real_)
    mean(logcpm274[idx274[i], gse274_meta$genotype == "P1-1"])
  }, numeric(1)),
  mean_logCPM_P1_2 = vapply(seq_along(candidate_genes), function(i) {
    if (is.na(idx274[i])) return(NA_real_)
    mean(logcpm274[idx274[i], gse274_meta$genotype == "P1-2"])
  }, numeric(1)),
  logFC_P1_2_minus_P1_1 = vapply(seq_along(candidate_genes), function(i) {
    if (is.na(idx274[i])) return(NA_real_)
    as.numeric(tab274$logFC[match(candidate_genes[i], tab274$gene)])
  }, numeric(1)),
  P_value = vapply(seq_along(candidate_genes), function(i) {
    if (is.na(idx274[i])) return(NA_real_)
    as.numeric(tab274$PValue[match(candidate_genes[i], tab274$gene)])
  }, numeric(1)),
  FDR = vapply(seq_along(candidate_genes), function(i) {
    if (is.na(idx274[i])) return(NA_real_)
    as.numeric(tab274$FDR[match(candidate_genes[i], tab274$gene)])
  }, numeric(1))
)
res274[, candidate_family_tested_n := sum(!is.na(P_value))]
res274[, FDR_candidate_tested := p.adjust(P_value, method = "BH")]
write_table(gse274_meta, file.path(out_dir, "GSE274269_sample_metadata.csv"))
write_table(res274, file.path(out_dir, "GSE274269_candidate_nodes.csv"))

## GSE179051: HeLa cells, M. pneumoniae (n=2) versus PBS (n=2).
gse179_counts <- fread(
  cmd = paste("gzip -cd", shQuote(file.path(base_dir, "GSE179051_mycoplasma_genecode_counts.txt.gz"))),
  skip = 1L,
  data.table = TRUE
)
stopifnot(identical(names(gse179_counts)[1:6], c("Geneid", "Chr", "Start", "End", "Strand", "Length")))
sample_names <- names(gse179_counts)[7:12]
valid_rows179 <- grepl("^ENSG[0-9]+$", gse179_counts$Geneid) &
  complete.cases(gse179_counts[, 7:12, with = FALSE])
gse179_counts <- gse179_counts[valid_rows179]
counts179 <- do.call(cbind, lapply(gse179_counts[, 7:12, with = FALSE], as.numeric))
colnames(counts179) <- sample_names
rownames(counts179) <- gse179_counts$Geneid
storage.mode(counts179) <- "numeric"

keys <- sub("\\..*$", "", rownames(counts179))
mapped <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(keys),
  keytype = "ENSEMBL",
  columns = "SYMBOL"
)
mapped <- mapped[!is.na(mapped$SYMBOL) & mapped$SYMBOL != "", c("ENSEMBL", "SYMBOL")]
mapped <- mapped[!duplicated(mapped$ENSEMBL), ]
symbol_for_row <- mapped$SYMBOL[match(keys, mapped$ENSEMBL)]
keep <- !is.na(symbol_for_row)
counts179 <- counts179[keep, , drop = FALSE]
rownames(counts179) <- symbol_for_row[keep]
if (anyDuplicated(rownames(counts179))) {
  count_dt <- as.data.table(counts179, keep.rownames = "gene")
  count_dt <- count_dt[, lapply(.SD, sum), by = gene]
  counts179 <- as.matrix(count_dt[, -1L, with = FALSE])
  rownames(counts179) <- count_dt$gene
  storage.mode(counts179) <- "numeric"
}

gse179_meta <- data.table(
  sample = sample_names,
  condition = factor(c("M_genitalium", "M_genitalium", "M_pneumoniae", "M_pneumoniae", "PBS", "PBS"),
                     levels = c("PBS", "M_pneumoniae", "M_genitalium"))
)
design179 <- model.matrix(~ 0 + condition, data = gse179_meta)
colnames(design179) <- levels(gse179_meta$condition)
y179_all <- DGEList(counts = counts179, samples = gse179_meta)
y179_all <- calcNormFactors(y179_all)
logcpm179 <- cpm(y179_all, log = TRUE, prior.count = 1)
keep179 <- filterByExpr(y179_all, design = design179)
y179 <- y179_all[keep179, , keep.lib.sizes = FALSE]
y179 <- calcNormFactors(y179)
fit179 <- glmQLFit(y179, design179, robust = TRUE)
contrast179 <- makeContrasts(M_pneumoniae - PBS, levels = design179)
qlf179 <- glmQLFTest(fit179, contrast = contrast179)
tab179 <- topTags(qlf179, n = Inf, sort.by = "none")$table
tab179$gene <- rownames(tab179)

idx179 <- match(candidate_genes, rownames(logcpm179))
res179 <- data.table(
  dataset = "GSE179051",
  comparison = "M_pneumoniae_vs_PBS_in_HeLa",
  gene = candidate_genes,
  detected_n = vapply(seq_along(candidate_genes), function(i) {
    if (is.na(idx179[i])) return(0L)
    sum(counts179[idx179[i], ] > 0)
  }, integer(1)),
  detected_rate = vapply(seq_along(candidate_genes), function(i) {
    if (is.na(idx179[i])) return(0)
    mean(counts179[idx179[i], ] > 0)
  }, numeric(1)),
  mean_logCPM_PBS = vapply(seq_along(candidate_genes), function(i) {
    if (is.na(idx179[i])) return(NA_real_)
    mean(logcpm179[idx179[i], gse179_meta$condition == "PBS"])
  }, numeric(1)),
  mean_logCPM_M_pneumoniae = vapply(seq_along(candidate_genes), function(i) {
    if (is.na(idx179[i])) return(NA_real_)
    mean(logcpm179[idx179[i], gse179_meta$condition == "M_pneumoniae"])
  }, numeric(1)),
  logFC_M_pneumoniae_minus_PBS = vapply(seq_along(candidate_genes), function(i) {
    if (is.na(idx179[i])) return(NA_real_)
    as.numeric(tab179$logFC[match(candidate_genes[i], tab179$gene)])
  }, numeric(1)),
  P_value = vapply(seq_along(candidate_genes), function(i) {
    if (is.na(idx179[i])) return(NA_real_)
    as.numeric(tab179$PValue[match(candidate_genes[i], tab179$gene)])
  }, numeric(1)),
  FDR = vapply(seq_along(candidate_genes), function(i) {
    if (is.na(idx179[i])) return(NA_real_)
    as.numeric(tab179$FDR[match(candidate_genes[i], tab179$gene)])
  }, numeric(1))
)
res179[, candidate_family_tested_n := sum(!is.na(P_value))]
res179[, FDR_candidate_tested := p.adjust(P_value, method = "BH")]
write_table(gse179_meta, file.path(out_dir, "GSE179051_sample_metadata.csv"))
write_table(res179, file.path(out_dir, "GSE179051_candidate_nodes.csv"))

summary <- data.table(
  dataset = c("GSE274269", "GSE179051"),
  samples = c(ncol(count_mat), ncol(counts179)),
  sample_context = c("M. pneumoniae-positive BALF total RNA; P1-1/P1-2 only", "HeLa; M. pneumoniae/PBS; n=2 per condition"),
  analysis_role = c("Independent BALF detectability and P1 genotype exploration", "Pathogen-stimulation direction check"),
  clinical_severity_available = c(FALSE, FALSE),
  cell_type_available = c(FALSE, FALSE)
)
write_table(summary, file.path(out_dir, "public_database_validation_summary.csv"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

cat("Completed public database small validation.\n")
cat("GSE274269 candidates:", nrow(res274), "\n")
cat("GSE179051 candidates:", nrow(res179), "\n")
