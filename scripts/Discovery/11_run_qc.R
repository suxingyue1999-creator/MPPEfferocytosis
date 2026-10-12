#!/usr/bin/env Rscript

local_lib <- Sys.getenv("BD_R_LIB", "inputs/discovery/bd_rhapsody_10/downstream/10samples_scrna/env/host-Rlib")
.libPaths(c(local_lib, .libPaths()))
suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(SoupX)
  library(DropletUtils)
  library(ggplot2)
})

set.seed(20260726)
samples <- c("D01", "D02", "D03", "D04", "D05", "D06", "D07", "D08", "D09", "D10")
input_root <- Sys.getenv("BD_INPUT_ROOT", "inputs/discovery/bd_rhapsody_10/work/outputs")
analysis_root <- Sys.getenv("BD_ANALYSIS_ROOT", "inputs/discovery/bd_rhapsody_10/downstream/10samples_scrna")
result_dir <- file.path(analysis_root, "results", "01_qc")
tmp_root <- file.path(analysis_root, "tmp", "01_qc")
state_dir <- file.path(analysis_root, "state")
emptydrops_dir <- file.path(result_dir, "emptydrops")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tmp_root, recursive = TRUE, showWarnings = FALSE)
dir.create(state_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(emptydrops_dir, recursive = TRUE, showWarnings = FALSE)

metadata <- read.csv(file.path(analysis_root, "results", "00_audit", "10samples_metadata_locked.csv"), stringsAsFactors = FALSE)
pipeline_metrics <- read.csv(file.path(analysis_root, "results", "00_audit", "pipeline_metrics_by_sample.csv"), stringsAsFactors = FALSE)
scrublet_summary <- read.csv(file.path(result_dir, "scrublet", "scrublet_summary.csv"), stringsAsFactors = FALSE)

read_mex <- function(sample_id, unfiltered = FALSE) {
  suffix <- if (unfiltered) "_RSEC_MolsPerCell_Unfiltered_MEX.zip" else "_RSEC_MolsPerCell_MEX.zip"
  role <- if (unfiltered) "unfiltered" else "filtered"
  archive <- file.path(input_root, sample_id, paste0(sample_id, suffix))
  target <- file.path(tmp_root, sample_id, role)
  unlink(target, recursive = TRUE, force = TRUE)
  dir.create(target, recursive = TRUE, showWarnings = FALSE)
  utils::unzip(archive, exdir = target)
  value <- Read10X(target, gene.column = 2, cell.column = 1, unique.features = TRUE)
  if (is.list(value)) value <- value[["Gene Expression"]]
  as(value, "dgCMatrix")
}

robust_log_bounds <- function(values, lower_floor, upper_mads = 4) {
  transformed <- log10(values + 1)
  center <- median(transformed)
  spread <- mad(transformed, constant = 1.4826)
  if (!is.finite(spread) || spread == 0) spread <- stats::sd(transformed)
  if (!is.finite(spread) || spread == 0) spread <- 0
  c(
    lower = max(lower_floor, floor(10^(center - 3 * spread) - 1)),
    upper = ceiling(10^(center + upper_mads * spread) - 1)
  )
}

safe_fraction <- function(counts, genes, totals) {
  present <- intersect(genes, rownames(counts))
  if (!length(present)) return(rep(0, ncol(counts)))
  as.numeric(Matrix::colSums(counts[present, , drop = FALSE])) / pmax(totals, 1) * 100
}

sample_workers <- suppressWarnings(as.integer(Sys.getenv("BD_SAMPLE_WORKERS", "2")))
if (!is.finite(sample_workers) || sample_workers < 1L) sample_workers <- 1L
sample_workers <- min(sample_workers, length(samples))

process_sample <- function(index) {
  sample_id <- samples[[index]]
  set.seed(20260726 + index)
  message("QC direct MEX: ", sample_id)
  counts <- read_mex(sample_id, unfiltered = FALSE)
  raw <- read_mex(sample_id, unfiltered = TRUE)
  barcodes <- colnames(counts)
  if (anyDuplicated(barcodes)) stop(sample_id, ": duplicated filtered barcodes")
  if (!all(barcodes %in% colnames(raw))) stop(sample_id, ": filtered barcodes absent from unfiltered MEX")
  if (!all(rownames(counts) %in% rownames(raw))) stop(sample_id, ": filtered features absent from unfiltered MEX")

  totals <- as.numeric(Matrix::colSums(counts))
  features <- as.numeric(Matrix::colSums(counts > 0))
  doublet_path <- file.path(result_dir, "scrublet", paste0(sample_id, "_scrublet.csv.gz"))
  doublets <- read.csv(gzfile(doublet_path), stringsAsFactors = FALSE)
  match_index <- match(barcodes, doublets$barcode)
  if (anyNA(match_index) || anyDuplicated(doublets$barcode)) stop(sample_id, ": Scrublet barcode mismatch")
  doublets <- doublets[match_index, , drop = FALSE]
  doublet_text <- tolower(trimws(as.character(doublets$predicted_doublet)))
  if (any(!doublet_text %in% c("true", "false"))) stop(sample_id, ": invalid Scrublet predicted_doublet values")
  doublets$predicted_doublet <- doublet_text == "true"

  percent_mt <- safe_fraction(counts, grep("^MT-", rownames(counts), value = TRUE, ignore.case = TRUE), totals)
  percent_ribo <- safe_fraction(counts, grep("^RP[SL]", rownames(counts), value = TRUE, ignore.case = TRUE), totals)
  percent_hb <- safe_fraction(counts, grep("^HB[ABDEGMQZ]", rownames(counts), value = TRUE, ignore.case = TRUE), totals)

  count_bounds <- robust_log_bounds(totals, lower_floor = 200)
  feature_bounds <- robust_log_bounds(features, lower_floor = 100)
  mt_center <- median(percent_mt)
  mt_spread <- mad(percent_mt, constant = 1.4826)
  if (!is.finite(mt_spread)) mt_spread <- 0
  mt_upper <- min(30, max(15, mt_center + 4 * mt_spread))

  message("emptyDrops: ", sample_id)
  ed <- emptyDrops(
    raw,
    lower = 100,
    niters = 10000,
    BPPARAM = BiocParallel::SerialParam()
  )
  ed_df <- as.data.frame(ed)
  ed_df$sample_id <- sample_id
  ed_df$barcode <- colnames(raw)
  ed_df$is_bd_filtered <- ed_df$barcode %in% barcodes
  ed_df$tested <- !is.na(ed_df$PValue)
  ed_df$emptydrops_fdr_0.01 <- ed_df$tested & !is.na(ed_df$FDR) & ed_df$FDR <= 0.01
  write.csv(
    ed_df[, c("sample_id", "barcode", "Total", "PValue", "FDR", "Limited", "tested", "emptydrops_fdr_0.01", "is_bd_filtered")],
    gzfile(file.path(emptydrops_dir, paste0(sample_id, "_emptydrops.csv.gz"))),
    row.names = FALSE
  )
  bd_tested <- ed_df$is_bd_filtered & ed_df$tested
  bd_confirmed <- ed_df$is_bd_filtered & ed_df$emptydrops_fdr_0.01
  emptydrops_only <- !ed_df$is_bd_filtered & ed_df$emptydrops_fdr_0.01
  ambient_background <- !ed_df$is_bd_filtered & ed_df$tested & !ed_df$emptydrops_fdr_0.01
  emptydrops_summary <- data.frame(
    sample_id = sample_id,
    unfiltered_barcodes = nrow(ed_df),
    tested_barcodes = sum(ed_df$tested),
    bd_filtered_cells = length(barcodes),
    bd_cells_tested = sum(bd_tested),
    bd_cells_confirmed_fdr_0.01 = sum(bd_confirmed),
    bd_cells_not_significant = sum(bd_tested & !ed_df$emptydrops_fdr_0.01),
    bd_cells_untested = sum(ed_df$is_bd_filtered & !ed_df$tested),
    emptydrops_only_candidates = sum(emptydrops_only),
    ambient_background_barcodes = sum(ambient_background),
    bd_tested_fraction = sum(bd_tested) / length(barcodes),
    bd_confirmed_fraction = sum(bd_confirmed) / length(barcodes),
    emptydrops_only_to_bd_ratio = sum(emptydrops_only) / length(barcodes),
    stringsAsFactors = FALSE
  )

  ranks <- barcodeRanks(raw, lower = 100)
  knee <- as.numeric(S4Vectors::metadata(ranks)$knee)
  inflection <- as.numeric(S4Vectors::metadata(ranks)$inflection)
  raw_for_soup <- raw[match(rownames(counts), rownames(raw)), , drop = FALSE]
  rownames(raw_for_soup) <- rownames(counts)
  toc <- raw_for_soup[, barcodes, drop = FALSE]

  prelim <- CreateSeuratObject(counts = counts, min.cells = 0, min.features = 0)
  prelim <- NormalizeData(prelim, verbose = FALSE)
  prelim <- FindVariableFeatures(prelim, selection.method = "vst", nfeatures = 1500, verbose = FALSE)
  prelim <- ScaleData(prelim, features = VariableFeatures(prelim), verbose = FALSE)
  prelim <- RunPCA(prelim, features = VariableFeatures(prelim), npcs = 20, approx = TRUE, verbose = FALSE)
  prelim <- FindNeighbors(prelim, dims = 1:15, verbose = FALSE)
  prelim <- FindClusters(prelim, resolution = 0.5, algorithm = 1, random.seed = 20260726 + index, verbose = FALSE)
  clusters <- as.character(Idents(prelim))
  names(clusters) <- colnames(prelim)

  soup <- SoupChannel(tod = raw_for_soup, toc = toc, calcSoupProfile = TRUE)
  soup <- setClusters(soup, clusters[colnames(toc)])
  soupx_error <- NA_character_
  soup <- tryCatch(
    autoEstCont(soup, doPlot = FALSE),
    error = function(e) {
      soupx_error <<- conditionMessage(e)
      soup
    }
  )
  soup_profile <- soup$soupProfile
  soup_profile$gene <- rownames(soup_profile)
  soup_profile <- soup_profile[order(soup_profile$est, decreasing = TRUE), , drop = FALSE]
  top_soup <- head(soup_profile, 20)
  ambient <- data.frame(
    sample_id = sample_id,
    rank = seq_len(nrow(top_soup)),
    gene = top_soup$gene,
    soup_fraction = top_soup$est,
    stringsAsFactors = FALSE
  )
  ambient_proxy <- safe_fraction(counts, top_soup$gene, totals)
  rho <- if ("rho" %in% colnames(soup$metaData)) median(soup$metaData$rho, na.rm = TRUE) else NA_real_
  if (!is.finite(rho)) rho <- NA_real_

  pass_low <- totals >= count_bounds[["lower"]] & features >= feature_bounds[["lower"]]
  pass_high <- totals <= count_bounds[["upper"]] & features <= feature_bounds[["upper"]]
  pass_mt <- percent_mt <= mt_upper
  pass_doublet <- !doublets$predicted_doublet
  proposed_keep <- pass_low & pass_high & pass_mt & pass_doublet

  qc <- data.frame(
    sample_id = sample_id,
    severity = metadata$severity_en[match(sample_id, metadata$sample_id)],
    barcode = barcodes,
    cell_id = paste(sample_id, barcodes, sep = "_"),
    nCount_RNA = totals,
    nFeature_RNA = features,
    percent_mt = percent_mt,
    percent_ribo = percent_ribo,
    percent_hb = percent_hb,
    doublet_score = doublets$doublet_score,
    predicted_doublet = doublets$predicted_doublet,
    ambient_top20_percent = ambient_proxy,
    emptydrops_tested = ed_df$tested[match(barcodes, ed_df$barcode)],
    emptydrops_fdr = ed_df$FDR[match(barcodes, ed_df$barcode)],
    emptydrops_significant = ed_df$emptydrops_fdr_0.01[match(barcodes, ed_df$barcode)],
    below_barcode_inflection = totals < inflection,
    pass_low_complexity = pass_low,
    pass_high_outlier = pass_high,
    pass_mito = pass_mt,
    pass_doublet = pass_doublet,
    qc_keep_proposed = proposed_keep,
    stringsAsFactors = FALSE
  )

  threshold <- data.frame(
    sample_id = sample_id,
    min_nCount_RNA = count_bounds[["lower"]],
    max_nCount_RNA = count_bounds[["upper"]],
    min_nFeature_RNA = feature_bounds[["lower"]],
    max_nFeature_RNA = feature_bounds[["upper"]],
    max_percent_mt = mt_upper,
    barcode_knee_umi = knee,
    barcode_inflection_umi = inflection,
    emptydrops_lower = 100,
    emptydrops_niters = 10000,
    emptydrops_fdr_diagnostic = 0.01,
    scrublet_threshold = scrublet_summary$scrublet_threshold[match(sample_id, scrublet_summary$sample_id)],
    stringsAsFactors = FALSE
  )

  summary <- data.frame(
    sample_id = sample_id,
    severity = metadata$severity_en[match(sample_id, metadata$sample_id)],
    cells_input = ncol(counts),
    cells_proposed_keep = sum(proposed_keep),
    proposed_retention_rate = mean(proposed_keep),
    median_nCount_RNA = median(totals),
    median_nFeature_RNA = median(features),
    median_percent_mt = median(percent_mt),
    median_ambient_top20_percent = median(ambient_proxy),
    predicted_doublets = sum(!pass_doublet),
    predicted_doublet_rate = mean(!pass_doublet),
    cells_below_barcode_inflection = sum(totals < inflection),
    soupx_rho = rho,
    soupx_status = ifelse(is.na(soupx_error), "PASS", "DIAGNOSTIC_ERROR"),
    soupx_detail = ifelse(is.na(soupx_error), "autoEstCont completed", soupx_error),
    stringsAsFactors = FALSE
  )

  sample_result <- list(
    qc = qc,
    summary = summary,
    threshold = threshold,
    emptydrops_summary = emptydrops_summary,
    ambient = ambient
  )
  rm(counts, raw, raw_for_soup, toc, prelim, soup, ranks, ed, ed_df)
  unlink(file.path(tmp_root, sample_id), recursive = TRUE, force = TRUE)
  gc(verbose = FALSE)
  sample_result
}

message("Stage 01 sample workers: ", sample_workers)
sample_results <- parallel::mclapply(
  seq_along(samples),
  process_sample,
  mc.cores = sample_workers,
  mc.preschedule = FALSE,
  mc.set.seed = TRUE
)
failed_workers <- vapply(sample_results, inherits, logical(1), what = "try-error")
if (any(failed_workers)) stop("Stage 01 sample worker failure: ", paste(samples[failed_workers], collapse = ","))
qc_rows <- lapply(sample_results, `[[`, "qc")
summary_rows <- lapply(sample_results, `[[`, "summary")
threshold_rows <- lapply(sample_results, `[[`, "threshold")
emptydrops_summary_rows <- lapply(sample_results, `[[`, "emptydrops_summary")
ambient_rows <- lapply(sample_results, `[[`, "ambient")

cell_qc <- do.call(rbind, qc_rows)
summary <- do.call(rbind, summary_rows)
thresholds <- do.call(rbind, threshold_rows)
ambient <- do.call(rbind, ambient_rows)
emptydrops_summary <- do.call(rbind, emptydrops_summary_rows)
summary <- merge(summary, emptydrops_summary, by = "sample_id", all.x = TRUE, sort = FALSE)
summary <- merge(summary, pipeline_metrics, by = "sample_id", all.x = TRUE, sort = FALSE)
summary <- summary[match(samples, summary$sample_id), , drop = FALSE]

outlier_columns <- c(
  "cells_input", "median_nCount_RNA", "median_nFeature_RNA", "median_percent_mt",
  "predicted_doublet_rate", "bd_confirmed_fraction", "pct_reads_useful", "sequencing_saturation_pct"
)
outlier_scores <- sapply(outlier_columns, function(column) {
  values <- summary[[column]]
  spread <- mad(values, constant = 1.4826, na.rm = TRUE)
  if (!is.finite(spread) || spread == 0) spread <- sd(values, na.rm = TRUE)
  if (!is.finite(spread) || spread == 0) return(rep(0, length(values)))
  (values - median(values, na.rm = TRUE)) / spread
})
summary$max_abs_robust_z <- apply(abs(outlier_scores), 1, function(value) max(value, na.rm = TRUE))
summary$sample_outlier_flag <- summary$max_abs_robust_z > 3.5

guard_rows <- list()
add_guard <- function(sample_id, guardrail, passed, observed, threshold, detail) {
  guard_rows[[length(guard_rows) + 1L]] <<- data.frame(
    sample_id = sample_id,
    guardrail = guardrail,
    status = ifelse(isTRUE(passed), "PASS", "FAIL"),
    observed = as.character(observed),
    threshold = threshold,
    detail = detail,
    stringsAsFactors = FALSE
  )
}
for (sample_id in samples) {
  row <- summary[summary$sample_id == sample_id, , drop = FALSE]
  add_guard(sample_id, "retention_rate", row$proposed_retention_rate >= 0.50, row$proposed_retention_rate, ">=0.50", "Proposed QC retention")
  add_guard(sample_id, "retained_cells", row$cells_proposed_keep >= 500, row$cells_proposed_keep, ">=500", "Minimum analyzable cells")
  add_guard(sample_id, "scrublet_doublet_rate", row$predicted_doublet_rate <= 0.15, row$predicted_doublet_rate, "<=0.15", "High doublet calls require review")
  add_guard(sample_id, "sample_robust_z", row$max_abs_robust_z <= 3.5, row$max_abs_robust_z, "<=3.5", "Cross-sample robust outlier score")
  add_guard(sample_id, "emptydrops_bd_tested", row$bd_tested_fraction >= 0.95, row$bd_tested_fraction, ">=0.95", "BD cells evaluable by emptyDrops")
  add_guard(sample_id, "emptydrops_bd_concordance", row$bd_confirmed_fraction >= 0.80, row$bd_confirmed_fraction, ">=0.80", "Diagnostic concordance only; BD calls remain primary")
  add_guard(sample_id, "emptydrops_only_ratio", row$emptydrops_only_to_bd_ratio <= 1.00, row$emptydrops_only_to_bd_ratio, "<=1.00", "Large emptyDrops-only population requires review")
  add_guard(sample_id, "soupx_completed", row$soupx_status == "PASS", row$soupx_status, "PASS", row$soupx_detail)
  add_guard(sample_id, "soupx_rho", is.finite(row$soupx_rho) && row$soupx_rho <= 0.35, row$soupx_rho, "<=0.35", "High ambient RNA estimate requires review")
  add_guard(sample_id, "median_percent_mt", row$median_percent_mt <= 20, row$median_percent_mt, "<=20", "High sample-level mitochondrial fraction")
}
guardrails <- do.call(rbind, guard_rows)
guardrails <- guardrails[order(match(guardrails$sample_id, samples), guardrails$guardrail), , drop = FALSE]
failures <- guardrails[guardrails$status == "FAIL", , drop = FALSE]

write.csv(summary, file.path(result_dir, "10samples_qc_summary_pre_filter.csv"), row.names = FALSE)
write.csv(thresholds, file.path(result_dir, "qc_thresholds_proposed.csv"), row.names = FALSE)
write.csv(guardrails, file.path(result_dir, "qc_guardrails.csv"), row.names = FALSE)
write.csv(emptydrops_summary, file.path(result_dir, "emptydrops_summary.csv"), row.names = FALSE)
write.csv(ambient, file.path(result_dir, "ambient_top_genes.csv"), row.names = FALSE)
write.csv(cell_qc, gzfile(file.path(result_dir, "cell_qc_pre_filter.csv.gz")), row.names = FALSE)

plot_data <- rbind(
  data.frame(sample_id = cell_qc$sample_id, severity = cell_qc$severity, metric = "log10 UMI", value = log10(cell_qc$nCount_RNA + 1)),
  data.frame(sample_id = cell_qc$sample_id, severity = cell_qc$severity, metric = "log10 genes", value = log10(cell_qc$nFeature_RNA + 1)),
  data.frame(sample_id = cell_qc$sample_id, severity = cell_qc$severity, metric = "mitochondrial %", value = cell_qc$percent_mt),
  data.frame(sample_id = cell_qc$sample_id, severity = cell_qc$severity, metric = "Scrublet score", value = cell_qc$doublet_score)
)
p1 <- ggplot(plot_data, aes(sample_id, value, fill = severity)) +
  geom_violin(scale = "width", trim = TRUE) +
  facet_wrap(~metric, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c(mild = "#2878B5", severe = "#D64B40")) +
  theme_bw(base_size = 11) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(result_dir, "qc_distributions.pdf"), p1, width = 13, height = 9)
ggsave(file.path(result_dir, "qc_distributions.png"), p1, width = 13, height = 9, dpi = 160)

p2 <- ggplot(cell_qc, aes(nCount_RNA, nFeature_RNA, color = percent_mt)) +
  geom_point(alpha = 0.25, size = 0.35) + scale_x_log10() + scale_y_log10() +
  scale_color_viridis_c(option = "C") + facet_wrap(~sample_id, ncol = 5) +
  theme_bw(base_size = 10)
ggsave(file.path(result_dir, "qc_count_feature_scatter.pdf"), p2, width = 15, height = 7)

p3 <- ggplot(summary, aes(sample_id, proposed_retention_rate, fill = severity)) +
  geom_col(width = 0.72) + geom_text(aes(label = cells_proposed_keep), vjust = -0.2, size = 3) +
  scale_fill_manual(values = c(mild = "#2878B5", severe = "#D64B40")) +
  scale_y_continuous(limits = c(0, 1.08)) + theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(result_dir, "qc_retention.pdf"), p3, width = 10, height = 6)

concordance_plot <- ggplot(emptydrops_summary, aes(sample_id, bd_confirmed_fraction)) +
  geom_col(fill = "#3A7D44", width = 0.72) + geom_hline(yintercept = 0.80, linetype = 2, color = "#B33A3A") +
  scale_y_continuous(limits = c(0, 1)) + theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(result_dir, "emptydrops_concordance.pdf"), concordance_plot, width = 10, height = 6)
ggsave(file.path(result_dir, "emptydrops_concordance.png"), concordance_plot, width = 10, height = 6, dpi = 160)

approval_mode <- if (nrow(failures) == 0) "AUTO_GUARDRAILS" else "PENDING_EXCEPTION_REVIEW"
approval <- data.frame(
  approval_mode = approval_mode,
  status = ifelse(nrow(failures) == 0, "APPROVED", "PENDING"),
  failed_guardrails = nrow(failures),
  timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  stringsAsFactors = FALSE
)
write.csv(approval, file.path(result_dir, "qc_approval.csv"), row.names = FALSE)

if (nrow(failures) == 0) {
  approved <- thresholds
  approved$approval_mode <- "AUTO_GUARDRAILS"
  approved$approved_at <- approval$timestamp
  write.csv(approved, file.path(result_dir, "qc_thresholds_approved.csv"), row.names = FALSE)
  write.csv(approved, file.path(result_dir, "qc_thresholds.csv"), row.names = FALSE)
  cell_qc$qc_keep <- cell_qc$qc_keep_proposed
  summary$cells_qc_keep <- summary$cells_proposed_keep
  summary$retention_rate <- summary$proposed_retention_rate
  write.csv(summary, file.path(result_dir, "10samples_qc_summary.csv"), row.names = FALSE)
  write.csv(cell_qc, gzfile(file.path(result_dir, "cell_qc.csv.gz")), row.names = FALSE)
  keep_columns <- cell_qc[cell_qc$qc_keep, c("sample_id", "barcode", "cell_id"), drop = FALSE]
  write.csv(keep_columns, gzfile(file.path(result_dir, "cells_keep.csv.gz")), row.names = FALSE)
} else {
  pending_path <- file.path(state_dir, "01_qc_review.PENDING")
  pending_lines <- c(
    paste0("pending=", approval$timestamp),
    paste0("failed_guardrails=", nrow(failures)),
    apply(failures, 1, function(row) paste(row[["sample_id"]], row[["guardrail"]], row[["observed"]], row[["threshold"]], sep = "\t"))
  )
  writeLines(pending_lines, pending_path)
}

review_report <- c(
  "# Stage 01 QC threshold review packet",
  "",
  paste0("- Approval mode: ", approval_mode),
  paste0("- Input BD putative cells: ", format(sum(summary$cells_input), big.mark = ",")),
  paste0("- Cells retained by proposed thresholds: ", format(sum(summary$cells_proposed_keep), big.mark = ",")),
  paste0("- Failed guardrails: ", nrow(failures)),
  paste0("- Parallelism: ", sample_workers, " independent sample workers; emptyDrops remains serial within each sample."),
  "- Canonical counts: filtered RSEC MEX, read directly for every sample.",
  "- emptyDrops: lower=100, niters=10000, FDR<=0.01; diagnostic only.",
  "- BD filtered calls remain primary; emptyDrops-only barcodes are not added automatically.",
  "- SoupX estimates are diagnostic; primary counts are not altered.",
  "- Severity is not used to derive QC thresholds.",
  "- Stage 02 starts only after automatic guardrail approval or an explicit reviewed approval.",
  ""
)
if (nrow(failures)) {
  review_report <- c(review_report, "## Exceptions", "", apply(failures, 1, function(row) {
    paste0("- ", row[["sample_id"]], ": ", row[["guardrail"]], " observed=", row[["observed"]], " threshold=", row[["threshold"]])
  }))
}
writeLines(review_report, file.path(result_dir, "qc_review_packet.md"))
writeLines(review_report, file.path(result_dir, "qc_report.md"))
writeLines(capture.output(sessionInfo()), file.path(result_dir, "r_sessionInfo.txt"))
