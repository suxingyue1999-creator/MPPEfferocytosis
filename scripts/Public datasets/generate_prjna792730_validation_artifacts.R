#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

base_dir <- Sys.getenv("MPP_PUBLIC_ROOT", "inputs/public")
result_dir <- file.path(base_dir, "results", "PRJNA792730")
analysis_dir <- file.path(result_dir, "analysis")
figure_dir <- file.path(base_dir, "figures")
table_dir <- file.path(result_dir, "tables")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

count_path <- file.path(result_dir, "PRJNA792730_salmon_gene_counts.tsv")
qc_path <- file.path(result_dir, "PRJNA792730_salmon_qc.tsv")
candidate_path <- file.path(analysis_dir, "PRJNA792730_candidate_nodes.csv")
all_result_path <- file.path(analysis_dir, "PRJNA792730_all_gene_results.csv")
concordance_path <- file.path(analysis_dir, "PRJNA792730_publication_concordance_audit.csv")
stopifnot(file.exists(count_path), file.exists(qc_path), file.exists(candidate_path),
          file.exists(all_result_path), file.exists(concordance_path))

counts_dt <- fread(count_path)
qc <- fread(qc_path)
candidate <- fread(candidate_path)
all_results <- fread(all_result_path)
concordance <- fread(concordance_path)
fastq_audit_path <- file.path(base_dir, "results", "PRJNA792730_fastq_full_audit.tsv")
stopifnot(file.exists(fastq_audit_path))
fastq_audit <- fread(fastq_audit_path)
expected_samples <- c("C1-0h", "C2-0h", "C3-0h", "T1-12h", "T2-12h", "T3-12h",
                      "T1-24h", "T2-24h", "T3-24h")
expected_genes <- c("HIF1A", "THBS1", "CD36", "TGFB1", "MAPK1", "MAPK3", "SMAD2", "SMAD4",
                    "MERTK", "AXL", "C1QA", "C1QB", "C1QC", "APOE", "TREM2", "CTSD",
                    "LAMP1", "NPC2", "LIPA")
stopifnot(identical(names(counts_dt)[3:11], expected_samples), nrow(qc) == 9L,
          all(qc$inferred_library_type == "ISR"), nrow(candidate) == 57L,
          all(candidate$gene %in% expected_genes), all(unique(candidate$comparison) %in%
            c("12h_vs_0h", "24h_vs_0h", "24h_vs_12h")))

qc_export <- copy(qc)
qc_export[, `:=`(reads_processed = num_processed, reads_mapped = num_mapped,
                 mapping_rate_percent = 100 * mapping_rate,
                 inferred_library_type = as.character(inferred_library_type),
                 timepoint = sub("h$", " h", timepoint))]
qc_export[, c("num_processed", "num_mapped", "mapping_rate") := NULL]
setcolorder(qc_export, c("run", "sample", "timepoint", "reads_processed", "reads_mapped",
                         "mapping_rate_percent", "inferred_library_type", "salmon_version",
                         "gencode_target_mapping_rate"))
fwrite(qc_export, file.path(table_dir, "S_PRJNA792730_sample_QC.csv"), na = "NA")

timepoint_summary <- qc[, .(
  n_samples = .N,
  reads_processed_total = sum(num_processed),
  reads_mapped_total = sum(num_mapped),
  mapping_rate_min_percent = 100 * min(mapping_rate),
  mapping_rate_max_percent = 100 * max(mapping_rate)
), by = .(timepoint)]
setorder(timepoint_summary, timepoint)
fwrite(timepoint_summary, file.path(table_dir, "S_PRJNA792730_timepoint_summary.csv"), na = "NA")

candidate_export <- copy(candidate)
candidate_export[, status := fifelse(is.na(P_value),
  "not tested after expression filtering", "tested in the candidate-node family")]
candidate_export[, source_note := "Salmon gene-level estimated counts; edgeR quasi-likelihood; TMM normalization"]
fwrite(candidate_export, file.path(table_dir, "S_PRJNA792730_candidate_nodes.csv"), na = "NA")
fwrite(all_results, file.path(table_dir, "S_PRJNA792730_all_gene_results.csv"), na = "NA")
fwrite(concordance, file.path(table_dir, "S_PRJNA792730_publication_concordance.csv"), na = "NA")
fwrite(fastq_audit, file.path(table_dir, "S_PRJNA792730_fastq_audit.csv"), na = "NA")

candidate_family_summary <- candidate[, .(
  prespecified_nodes = 19L,
  tested_nodes = sum(!is.na(P_value)),
  not_tested_after_expression_filter = sum(is.na(P_value)),
  candidate_FDR_scope = "BH adjustment among tested candidate nodes in this comparison"
), by = .(comparison)]
fwrite(candidate_family_summary,
       file.path(table_dir, "S_PRJNA792730_candidate_family_summary.csv"), na = "NA")

count_mat <- as.matrix(counts_dt[, expected_samples, with = FALSE])
rownames(count_mat) <- counts_dt$gene_id
storage.mode(count_mat) <- "numeric"
sample_info <- data.table(
  sample = expected_samples,
  timepoint = factor(c("0 h", "0 h", "0 h", "12 h", "12 h", "12 h", "24 h", "24 h", "24 h"),
                     levels = c("0 h", "12 h", "24 h"))
)
y <- calcNormFactors(DGEList(counts = count_mat), method = "TMM")
logcpm <- cpm(y, log = TRUE, prior.count = 1)
gene_id_by_name <- counts_dt[!duplicated(gene_name), setNames(gene_id, gene_name)]
candidate_ids <- unname(gene_id_by_name[expected_genes])
candidate_ids <- candidate_ids[!is.na(candidate_ids)]
expr <- logcpm[intersect(candidate_ids, rownames(logcpm)), , drop = FALSE]
expr_names <- names(gene_id_by_name)[match(rownames(expr), gene_id_by_name)]
rownames(expr) <- expr_names
expr <- expr[rev(expected_genes[expected_genes %in% rownames(expr)]), , drop = FALSE]
expr_long <- as.data.table(as.table(expr))
setnames(expr_long, c("gene", "sample", "logCPM"))
expr_long[, gene := factor(as.character(gene), levels = rev(expected_genes))]
expr_long[, sample := factor(as.character(sample), levels = expected_samples)]

set3 <- c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462",
          "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")
time_cols <- c("0 h" = set3[1], "12 h" = set3[3], "24 h" = set3[5])

p_qc <- ggplot(qc_export, aes(x = sample, y = mapping_rate_percent, fill = timepoint)) +
  geom_col(width = 0.68, colour = "black", linewidth = 0.18) +
  geom_text(aes(label = sprintf("%.1f", mapping_rate_percent)), vjust = -0.35,
            family = "Times New Roman", size = 2.8) +
  scale_fill_manual(values = time_cols, drop = FALSE) +
  scale_y_continuous(limits = c(0, 50), breaks = c(0, 10, 20, 30, 40, 50),
                     labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0, 0.08))) +
  labs(x = NULL, y = "Salmon mapping rate") +
  theme_classic(base_family = "Times New Roman", base_size = 10) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
        axis.text.y = element_text(colour = "black"), axis.title = element_text(colour = "black"),
        plot.margin = margin(8, 10, 6, 8))

p_pca <- {
  pca <- prcomp(t(logcpm), scale. = TRUE)
  pca_dt <- data.table(sample = rownames(pca$x), PC1 = pca$x[, 1], PC2 = pca$x[, 2])
  pca_dt <- merge(pca_dt, sample_info, by = "sample", sort = FALSE)
  var_explained <- 100 * summary(pca)$importance[2, 1:2]
  ggplot(pca_dt, aes(PC1, PC2, fill = timepoint)) +
    geom_point(shape = 21, size = 3.2, colour = "black", stroke = 0.25) +
    scale_fill_manual(values = time_cols, drop = FALSE) +
    labs(x = sprintf("PC1 (%.1f%%)", var_explained[1]), y = sprintf("PC2 (%.1f%%)", var_explained[2])) +
    theme_classic(base_family = "Times New Roman", base_size = 10) +
  theme(legend.position = "bottom", legend.title = element_blank(),
          axis.text = element_text(colour = "black"), axis.title = element_text(colour = "black"),
          plot.margin = margin(8, 10, 6, 8))
}

p_heat <- ggplot(expr_long, aes(x = sample, y = gene, fill = logCPM)) +
  geom_tile(colour = "white", linewidth = 0.18) +
  scale_fill_gradientn(colours = c(set3[4], "white", set3[1], set3[7]), name = "logCPM") +
  labs(x = NULL, y = NULL) +
  theme_classic(base_family = "Times New Roman", base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", size = 7.5),
        axis.text.y = element_text(colour = "black", size = 7.5), axis.title = element_text(colour = "black"),
        legend.position = "right", plot.margin = margin(8, 10, 6, 8))

primary <- candidate[comparison == "24h_vs_0h"]
primary[, gene := factor(gene, levels = rev(expected_genes))]
p_primary <- ggplot(primary, aes(x = logFC, y = gene, colour = FDR_candidate_tested < 0.05)) +
  geom_vline(xintercept = 0, colour = "black", linewidth = 0.3) +
  geom_segment(aes(x = 0, xend = logFC, y = gene, yend = gene),
               linewidth = 0.45, colour = "#6B6B6B", na.rm = TRUE) +
  geom_point(size = 2.7, na.rm = TRUE) +
  scale_colour_manual(values = c(`TRUE` = set3[4], `FALSE` = set3[5]),
                      labels = c(`TRUE` = "candidate BH-FDR < 0.05", `FALSE` = "candidate BH-FDR >= 0.05"),
                      na.value = set3[9], name = NULL) +
  labs(x = "log2 fold change: 24 h minus 0 h", y = NULL) +
  theme_classic(base_family = "Times New Roman", base_size = 9.5) +
  theme(axis.text.y = element_text(colour = "black", size = 7.5), axis.text.x = element_text(colour = "black"),
        axis.title = element_text(colour = "black"), legend.position = "bottom",
        legend.text = element_text(size = 8), plot.margin = margin(10, 12, 8, 28))

figure <- (p_qc | p_pca) / (p_heat | p_primary) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(family = "Times New Roman", face = "bold", size = 12,
                                colour = "black"), plot.tag.position = c(0.01, 0.99))

stem <- file.path(figure_dir, "Supplementary_Figure_PRJNA792730_time_response")
ggsave(paste0(stem, ".svg"), figure, width = 11.6, height = 8.5, units = "in", device = "svg", bg = "white")
system2(file.path(base_dir, "render_prjna792730_pdf.sh"))
ggsave(paste0(stem, ".png"), figure, width = 11.6, height = 8.5, units = "in", dpi = 300, bg = "white")
ggsave(paste0(stem, ".tiff"), figure, width = 11.6, height = 8.5, units = "in", dpi = 600,
       compression = "lzw", bg = "white")
writeLines(capture.output(sessionInfo()), file.path(figure_dir, "Supplementary_Figure_PRJNA792730_time_response_sessionInfo.txt"))

primary <- candidate[comparison == "24h_vs_0h"]
primary_tested <- sum(!is.na(primary$P_value))
qc_min <- min(qc$mapping_rate)
qc_max <- max(qc$mapping_rate)
audit_lines <- c(
  "# PRJNA792730 external validation audit",
  "",
  "## Data identity",
  "",
  "The dataset contains nine paired-end bulk RNA-seq libraries from THP-1-derived macrophage-like cells infected with Mycoplasma pneumoniae: three libraries at 0 h, three at 12 h, and three at 24 h. It is an in-vitro time-response resource and is not a pediatric, BALF, clinical, severity, or patient-paired validation cohort.",
  "",
  "## Read and quantification audit",
  "",
  sprintf("All 18 FASTQ files passed the completed audit for expected byte size, ENA MD5, gzip integrity, FASTQ record counts, malformed records, and paired-read identifiers. Salmon 1.10.3 with GENCODE v44 produced complete output for all nine samples. The inferred library type was ISR in all nine samples. Mapping rates ranged from %.3f to %.3f (%.1f%% to %.1f%%).", qc_min, qc_max, 100 * qc_min, 100 * qc_max),
  "",
  "## Statistical audit",
  "",
  sprintf("Gene-level estimated counts were summarized by gene_id and analyzed with edgeR quasi-likelihood modeling, TMM normalization, and a no-intercept timepoint design. The three prespecified comparisons were 12 h versus 0 h, 24 h versus 0 h, and 24 h versus 12 h. Each comparison tested %d genes after filtering.", unique(all_results$tested_genes[all_results$comparison == "24h_vs_0h"])),
  sprintf("The 24 h versus 0 h comparison tested %d of the 19 prespecified candidate nodes after expression filtering. C1QA, C1QB, and C1QC were not tested in this model because they did not pass the expression filter; their missing statistics represent unavailable tests, not null effects.", primary_tested),
  "Candidate-node P values were adjusted separately within the prespecified 19-node family, with the effective tested family size reported in the candidate table. Whole-transcriptome FDR and candidate-family BH-FDR were kept as separate quantities.",
  "",
  "## Comparison with the source publication",
  "",
  sprintf("The source publication reported 927 upregulated and 1,162 downregulated genes (2,089 total) for its 24 h infection comparison. The present reanalysis reported %d genes with absolute log2 fold change > 1 and whole-transcriptome FDR < 0.05 (%d upregulated and %d downregulated). These totals are retained as a method-specific concordance audit and are not merged with the source publication result.", concordance$reanalysis_total, concordance$reanalysis_up, concordance$reanalysis_down),
  "",
  "## Interpretation boundary",
  "",
  "The public dataset can provide independent support for a time-dependent response of selected candidate nodes in an infected macrophage-like cell model. It cannot establish a pediatric severe-disease phenotype, a BALF cellular state, a clinical outcome association, a patient-level cross-layer relationship, or MEK/ERK phosphorylation. MAPK1 and MAPK3 remain transcript-level candidate relays only."
)
writeLines(audit_lines, file.path(result_dir, "PRJNA792730_validation_audit.md"))

primary_concordance <- concordance[comparison == "24h_vs_0h"]
audit_checks <- data.table(
  check_id = c("FASTQ_FILES", "FASTQ_RUNS", "FASTQ_INTEGRITY", "SALMON_OUTPUTS",
               "LIBRARY_TYPE", "MAPPING_RATE", "GENE_TESTS", "CANDIDATE_NODES",
               "CANDIDATE_UNTESTED", "PRIMARY_REANALYSIS", "SOURCE_PUBLICATION",
               "SAMPLE_IDENTITY"),
  domain = c("raw data", "raw data", "raw data", "quantification", "quantification",
             "quantification", "statistics", "statistics", "statistics", "concordance",
             "concordance", "interpretation"),
  observed = c(
    sprintf("%d of 18 FASTQ files present", nrow(fastq_audit)),
    sprintf("%d runs with both reads", uniqueN(fastq_audit$run)),
    sprintf("%d of %d files PASS for bytes, MD5, gzip, records, malformed records and paired IDs",
            sum(fastq_audit$status == "PASS"), nrow(fastq_audit)),
    "9 of 9 samples have complete Salmon output",
    sprintf("ISR in %d of 9 samples", sum(qc$inferred_library_type == "ISR")),
    sprintf("%.3f%% to %.3f%%", 100 * qc_min, 100 * qc_max),
    sprintf("19,518 genes in each of 3 prespecified comparisons"),
    "19 prespecified nodes; 16 tested in each comparison",
    "C1QA, C1QB and C1QC; NA denotes failure of expression filtering",
    sprintf("%d genes with |log2FC| > 1 and whole-transcriptome FDR < 0.05 (%d up, %d down)",
            primary_concordance$reanalysis_total, primary_concordance$reanalysis_up,
            primary_concordance$reanalysis_down),
    sprintf("%d genes reported in the source publication (%d up, %d down)",
            primary_concordance$publication_reported_total,
            primary_concordance$publication_reported_up,
            primary_concordance$publication_reported_down),
    "9 THP-1-derived macrophage-like libraries; 3 samples at each time point"
  ),
  expected = c("18", "9", "all PASS", "9", "all ISR", "reported range", "19,518 each",
               "19 prespecified", "documented NA rule", "method-specific result",
               "source-publication result", "time-response model"),
  status = c("PASS", "PASS", "PASS", "PASS", "PASS", "PASS", "PASS", "PASS",
             "PASS", "PASS", "PASS", "BOUNDARY"),
  source = c("PRJNA792730_fastq_full_audit.tsv", "PRJNA792730_fastq_full_audit.tsv",
             "PRJNA792730_fastq_full_audit.tsv", "PRJNA792730_salmon_qc.tsv",
             "PRJNA792730_salmon_qc.tsv", "PRJNA792730_salmon_qc.tsv",
             "PRJNA792730_all_gene_results.csv", "PRJNA792730_candidate_nodes.csv",
             "PRJNA792730_candidate_nodes.csv", "PRJNA792730_publication_concordance_audit.csv",
             "PRJNA792730_publication_concordance_audit.csv", "PRJNA792730_sample_metadata.csv")
)
fwrite(audit_checks, file.path(table_dir, "S_PRJNA792730_audit_checks.csv"), na = "NA")
cat("Generated PRJNA792730 validation tables, figure files, and audit report.\n")
