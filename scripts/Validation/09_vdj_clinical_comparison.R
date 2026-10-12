#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

set.seed(20260721)
options(stringsAsFactors = FALSE)

workspace <- Sys.getenv("MPP_WORKSPACE", "inputs/workspace")
work_root <- file.path(workspace, "工作文件", "23例单细胞重新分析")
input_dir <- file.path(work_root, "results", "04_vdj")
out_dir <- file.path(work_root, "results", "09_vdj_clinical_comparison")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

clinical_path <- file.path(workspace, "工作文件", "单细胞23例样本临床信息.csv")
metric_path <- file.path(input_dir, "sample_receptor_summary.csv")
cell_path <- file.path(input_dir, "cell_receptor_summary.tsv.gz")
contig_path <- file.path(input_dir, "vdj_contigs_selected.tsv.gz")
linked_path <- file.path(input_dir, "vdj_cell_level_with_metadata.tsv.gz")

required_paths <- c(clinical_path, metric_path, cell_path, contig_path, linked_path)
if (any(!file.exists(required_paths))) {
  stop("Missing required VDJ inputs: ", paste(required_paths[!file.exists(required_paths)], collapse = ", "))
}

read_tsv_gz <- function(path) {
  read.delim(
    gzfile(path), sep = "\t", quote = "", comment.char = "",
    check.names = FALSE, stringsAsFactors = FALSE
  )
}

write_tsv_gz <- function(data, path) {
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  write.table(data, con, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
}

adjust_na <- function(p) {
  out <- rep(NA_real_, length(p))
  keep <- is.finite(p)
  out[keep] <- p.adjust(p[keep], method = "BH")
  out
}

add_grouped_fdr <- function(data, group_cols, p_col = "adjusted_p_value", out_col = "FDR_BH") {
  if (!nrow(data)) {
    data[[out_col]] <- numeric()
    return(data)
  }
  keys <- interaction(data[, group_cols, drop = FALSE], drop = TRUE, lex.order = TRUE)
  data[[out_col]] <- ave(data[[p_col]], keys, FUN = adjust_na)
  data
}

gini_coefficient <- function(x) {
  x <- sort(as.numeric(x[x > 0]))
  n <- length(x)
  if (!n || sum(x) == 0) return(NA_real_)
  sum((2 * seq_len(n) - n - 1) * x) / (n * sum(x))
}

clamp01 <- function(x, eps = 1e-6) pmin(pmax(x, eps), 1 - eps)

transform_metric <- function(x, transformation) {
  if (transformation == "log1p") return(log1p(x))
  if (transformation == "log") return(log(pmax(x, 1e-8)))
  if (transformation == "logit") return(qlogis(clamp01(x)))
  x
}

clinical <- read.csv(clinical_path, check.names = FALSE, stringsAsFactors = FALSE)
clinical$sample_id <- clinical$sample_code
clinical$batch_wave <- as.character(clinical$batch_wave)
clinical <- unique(clinical[, c(
  "sample_id", "severity", "allergy_group", "plastic_group",
  "batch_wave"
)])
clinical$batch_wave <- factor(clinical$batch_wave, levels = c("2023_Dec", "2024_Spring"))
clinical$severity <- factor(clinical$severity, levels = c("轻症", "重症"))
clinical$allergy_group <- factor(clinical$allergy_group, levels = c("否", "是"))
clinical$plastic_group <- factor(clinical$plastic_group, levels = c("否", "是"))
clinical <- clinical[order(clinical$sample_id), ]

if (nrow(clinical) != 23L || anyDuplicated(clinical$sample_id)) {
  stop("Clinical design must contain exactly 23 unique samples")
}

design <- model.matrix(~ batch_wave + severity + allergy_group + plastic_group, data = clinical)
design_rank <- qr(design)$rank
if (design_rank != ncol(design)) stop("Clinical design matrix is not full rank")

design_export <- data.frame(
  clinical[, c("sample_id")],
  batch_wave = as.character(clinical$batch_wave),
  severity = as.character(clinical$severity),
  allergy_group = as.character(clinical$allergy_group),
  plastic_group = as.character(clinical$plastic_group),
  design,
  check.names = FALSE
)
write.csv(design_export, file.path(out_dir, "clinical_design_matrix.csv"), row.names = FALSE, fileEncoding = "UTF-8")

contrast_info <- list(
  severity = list(variable = "severity", coefficient = "severity重症", reference = "轻症", case = "重症", label = "重症_vs_轻症"),
  allergy = list(variable = "allergy_group", coefficient = "allergy_group是", reference = "否", case = "是", label = "过敏_vs_不过敏"),
  plastic = list(variable = "plastic_group", coefficient = "plastic_group是", reference = "否", case = "是", label = "塑形_vs_非塑形")
)

fit_continuous_comparisons <- function(data, value_col, strata_cols, transformation = "identity") {
  groups <- split(data, interaction(data[, strata_cols, drop = FALSE], drop = TRUE, lex.order = TRUE))
  rows <- list()
  index <- 0L
  for (current in groups) {
    current <- merge(
      current,
      clinical[, c("sample_id", "batch_wave", "severity", "allergy_group", "plastic_group")],
      by = "sample_id", all.x = TRUE, sort = FALSE
    )
    current <- current[is.finite(current[[value_col]]), , drop = FALSE]
    if (nrow(current) < 8L) next
    current$transformed_value <- transform_metric(current[[value_col]], transformation)
    fit <- tryCatch(
      lm(transformed_value ~ batch_wave + severity + allergy_group + plastic_group, data = current),
      error = function(e) NULL
    )
    if (is.null(fit) || fit$rank != length(coef(fit))) next
    coefficient_table <- summary(fit)$coefficients
    for (contrast_name in names(contrast_info)) {
      info <- contrast_info[[contrast_name]]
      coefficient <- info$coefficient
      if (!coefficient %in% rownames(coefficient_table)) next
      raw_groups <- as.character(current[[info$variable]])
      reference_values <- current[[value_col]][raw_groups == info$reference]
      case_values <- current[[value_col]][raw_groups == info$case]
      raw_p <- if (length(reference_values) && length(case_values)) {
        suppressWarnings(wilcox.test(case_values, reference_values, exact = FALSE)$p.value)
      } else NA_real_
      index <- index + 1L
      row <- data.frame(
        contrast = info$label,
        reference = info$reference,
        case = info$case,
        samples_used = nrow(current),
        reference_n = length(reference_values),
        case_n = length(case_values),
        reference_median = median(reference_values, na.rm = TRUE),
        case_median = median(case_values, na.rm = TRUE),
        raw_median_difference = median(case_values, na.rm = TRUE) - median(reference_values, na.rm = TRUE),
        transformation = transformation,
        adjusted_beta = coefficient_table[coefficient, "Estimate"],
        adjusted_standard_error = coefficient_table[coefficient, "Std. Error"],
        adjusted_p_value = coefficient_table[coefficient, "Pr(>|t|)"],
        unadjusted_wilcoxon_p_value = raw_p,
        stringsAsFactors = FALSE
      )
      strata_values <- current[1, strata_cols, drop = FALSE]
      rows[[index]] <- cbind(strata_values, row)
    }
  }
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

cell_receptor <- read_tsv_gz(cell_path)
cell_receptor$clone_size <- as.numeric(cell_receptor$clone_size)
cell_receptor$paired_chain <- as.logical(cell_receptor$paired_chain)

metric_rows <- lapply(split(
  cell_receptor,
  interaction(cell_receptor$sample_id, cell_receptor$receptor, drop = TRUE, lex.order = TRUE)
), function(df) {
  clone_counts <- as.numeric(table(df$clone_id_within_sample))
  clone_counts <- sort(clone_counts, decreasing = TRUE)
  total <- sum(clone_counts)
  probabilities <- clone_counts / total
  shannon <- -sum(probabilities * log(probabilities))
  simpson_diversity <- 1 - sum(probabilities^2)
  richness <- length(clone_counts)
  clones_to_50 <- which(cumsum(probabilities) >= 0.5)[1]
  data.frame(
    sample_id = unique(df$sample_id),
    receptor = unique(df$receptor),
    cells_with_productive_receptor = nrow(df),
    unique_clones = richness,
    paired_fraction = mean(df$paired_chain, na.rm = TRUE),
    singleton_fraction = mean(df$clone_size == 1),
    expanded_fraction_ge2 = mean(df$clone_size >= 2),
    expanded_fraction_ge5 = mean(df$clone_size >= 5),
    max_clone_size = max(clone_counts),
    max_clone_fraction = max(probabilities),
    top10_clone_fraction = sum(head(probabilities, 10)),
    mean_clone_size = mean(clone_counts),
    shannon = shannon,
    simpson_diversity = simpson_diversity,
    inv_simpson = 1 / sum(probabilities^2),
    pielou_evenness = if (richness > 1) shannon / log(richness) else 1,
    clonality = if (richness > 1) 1 - shannon / log(richness) else 0,
    gini = gini_coefficient(clone_counts),
    d50_clone_count = clones_to_50,
    d50_clone_fraction = clones_to_50 / richness,
    stringsAsFactors = FALSE
  )
})
metrics <- do.call(rbind, metric_rows)
metrics <- merge(metrics, clinical, by = "sample_id", all.x = TRUE, sort = FALSE)
metrics <- metrics[order(metrics$sample_id, metrics$receptor), ]
write.csv(metrics, file.path(out_dir, "repertoire_metrics_by_sample.csv"), row.names = FALSE, fileEncoding = "UTF-8")

metric_transformations <- c(
  cells_with_productive_receptor = "log1p",
  unique_clones = "log1p",
  paired_fraction = "logit",
  singleton_fraction = "logit",
  expanded_fraction_ge2 = "logit",
  expanded_fraction_ge5 = "logit",
  max_clone_size = "log1p",
  max_clone_fraction = "logit",
  top10_clone_fraction = "logit",
  mean_clone_size = "log",
  shannon = "identity",
  simpson_diversity = "logit",
  inv_simpson = "log",
  pielou_evenness = "logit",
  clonality = "logit",
  gini = "logit",
  d50_clone_count = "log1p",
  d50_clone_fraction = "logit"
)

metric_comparison_rows <- lapply(names(metric_transformations), function(metric) {
  current <- metrics[, c("sample_id", "receptor", metric), drop = FALSE]
  names(current)[3] <- "value"
  result <- fit_continuous_comparisons(
    current, value_col = "value", strata_cols = c("receptor"),
    transformation = metric_transformations[[metric]]
  )
  if (nrow(result)) result$metric <- metric
  result
})
metric_comparisons <- do.call(rbind, metric_comparison_rows)
metric_comparisons <- metric_comparisons[, c(
  "receptor", "metric", setdiff(names(metric_comparisons), c("receptor", "metric"))
)]
metric_comparisons <- add_grouped_fdr(metric_comparisons, "contrast")
metric_comparisons$unadjusted_wilcoxon_FDR_BH <- ave(
  metric_comparisons$unadjusted_wilcoxon_p_value,
  metric_comparisons$contrast,
  FUN = adjust_na
)
metric_comparisons <- metric_comparisons[order(metric_comparisons$contrast, metric_comparisons$FDR_BH), ]
write.csv(metric_comparisons, file.path(out_dir, "repertoire_metric_clinical_comparisons.csv"), row.names = FALSE, fileEncoding = "UTF-8")

contigs <- read_tsv_gz(contig_path)
contigs$productive <- as.logical(contigs$productive)
contigs$cdr3_aa <- as.character(contigs$cdr3_aa)
productive <- contigs[contigs$productive & !is.na(contigs$cdr3_aa) & nzchar(contigs$cdr3_aa), ]
productive$cdr3_aa_length <- as.numeric(productive$cdr3_aa_length)

cdr3_groups <- split(
  productive,
  interaction(productive$sample_id, productive$receptor, productive$locus, drop = TRUE, lex.order = TRUE)
)
cdr3_by_sample <- do.call(rbind, lapply(cdr3_groups, function(df) {
  values <- df$cdr3_aa_length[is.finite(df$cdr3_aa_length)]
  data.frame(
    sample_id = unique(df$sample_id),
    receptor = unique(df$receptor),
    locus = unique(df$locus),
    productive_contigs = length(values),
    mean_cdr3_aa_length = mean(values),
    median_cdr3_aa_length = median(values),
    p90_cdr3_aa_length = unname(quantile(values, 0.9)),
    stringsAsFactors = FALSE
  )
}))
cdr3_by_sample <- merge(cdr3_by_sample, clinical, by = "sample_id", all.x = TRUE, sort = FALSE)
cdr3_by_sample <- cdr3_by_sample[order(cdr3_by_sample$sample_id, cdr3_by_sample$receptor, cdr3_by_sample$locus), ]
write.csv(cdr3_by_sample, file.path(out_dir, "cdr3_length_metrics_by_sample.csv"), row.names = FALSE, fileEncoding = "UTF-8")

cdr3_metrics <- c("mean_cdr3_aa_length", "median_cdr3_aa_length", "p90_cdr3_aa_length")
cdr3_comparison_rows <- lapply(cdr3_metrics, function(metric) {
  current <- cdr3_by_sample[, c("sample_id", "receptor", "locus", metric), drop = FALSE]
  names(current)[4] <- "value"
  result <- fit_continuous_comparisons(
    current, value_col = "value", strata_cols = c("receptor", "locus"),
    transformation = "identity"
  )
  if (nrow(result)) result$metric <- metric
  result
})
cdr3_comparisons <- do.call(rbind, cdr3_comparison_rows)
cdr3_comparisons <- cdr3_comparisons[, c(
  "receptor", "locus", "metric",
  setdiff(names(cdr3_comparisons), c("receptor", "locus", "metric"))
)]
cdr3_comparisons <- add_grouped_fdr(cdr3_comparisons, "contrast")
cdr3_comparisons <- cdr3_comparisons[order(cdr3_comparisons$contrast, cdr3_comparisons$FDR_BH), ]
write.csv(cdr3_comparisons, file.path(out_dir, "cdr3_length_clinical_comparisons.csv"), row.names = FALSE, fileEncoding = "UTF-8")

analyse_segment_usage <- function(data, feature_cols, feature_type, min_total, min_samples) {
  data <- data[data$productive & data$receptor %in% c("BCR", "TCR"), , drop = FALSE]
  keep_complete <- Reduce(`&`, lapply(feature_cols, function(column) !is.na(data[[column]]) & nzchar(data[[column]])))
  data <- data[keep_complete, , drop = FALSE]
  if (!nrow(data)) return(list(comparisons = data.frame(), features = data.frame()))

  data$feature <- do.call(paste, c(data[feature_cols], sep = "|"))
  count_table <- aggregate(
    rep(1L, nrow(data)),
    by = list(
      sample_id = data$sample_id,
      receptor = data$receptor,
      locus = data$locus,
      feature = data$feature
    ),
    FUN = sum
  )
  names(count_table)[5] <- "feature_count"
  denominators <- aggregate(
    feature_count ~ sample_id + receptor + locus,
    data = count_table,
    FUN = sum
  )
  names(denominators)[4] <- "locus_total"

  feature_summary <- aggregate(
    feature_count ~ receptor + locus + feature,
    data = count_table,
    FUN = sum
  )
  names(feature_summary)[4] <- "total_count"
  prevalence <- aggregate(
    sample_id ~ receptor + locus + feature,
    data = count_table[count_table$feature_count > 0, ],
    FUN = function(x) length(unique(x))
  )
  names(prevalence)[4] <- "samples_present"
  feature_summary <- merge(feature_summary, prevalence, by = c("receptor", "locus", "feature"), all.x = TRUE)
  feature_summary$feature_type <- feature_type
  feature_summary$min_total_threshold <- min_total
  feature_summary$min_samples_threshold <- min_samples
  feature_summary$passed_filter <- feature_summary$total_count >= min_total & feature_summary$samples_present >= min_samples
  tested_features <- feature_summary[feature_summary$passed_filter, , drop = FALSE]

  rows <- list()
  index <- 0L
  for (i in seq_len(nrow(tested_features))) {
    feature_row <- tested_features[i, ]
    current_denominators <- denominators[
      denominators$receptor == feature_row$receptor & denominators$locus == feature_row$locus,
      , drop = FALSE
    ]
    current_counts <- count_table[
      count_table$receptor == feature_row$receptor &
        count_table$locus == feature_row$locus &
        count_table$feature == feature_row$feature,
      c("sample_id", "feature_count"), drop = FALSE
    ]
    current <- merge(current_denominators, current_counts, by = "sample_id", all.x = TRUE, sort = FALSE)
    current$feature_count[is.na(current$feature_count)] <- 0L
    current <- merge(
      current,
      clinical[, c("sample_id", "batch_wave", "severity", "allergy_group", "plastic_group")],
      by = "sample_id", all.x = TRUE, sort = FALSE
    )
    current <- current[current$locus_total > 0, , drop = FALSE]
    if (nrow(current) < 8L) next
    fit <- tryCatch(
      suppressWarnings(glm(
        cbind(feature_count, locus_total - feature_count) ~
          batch_wave + severity + allergy_group + plastic_group,
        family = quasibinomial(), data = current
      )),
      error = function(e) NULL
    )
    if (is.null(fit) || fit$rank != length(coef(fit))) next
    coefficient_table <- summary(fit)$coefficients
    current$feature_fraction <- current$feature_count / current$locus_total
    for (contrast_name in names(contrast_info)) {
      info <- contrast_info[[contrast_name]]
      coefficient <- info$coefficient
      if (!coefficient %in% rownames(coefficient_table)) next
      group_values <- as.character(current[[info$variable]])
      reference_fraction <- current$feature_fraction[group_values == info$reference]
      case_fraction <- current$feature_fraction[group_values == info$case]
      index <- index + 1L
      rows[[index]] <- data.frame(
        feature_type = feature_type,
        receptor = feature_row$receptor,
        locus = feature_row$locus,
        feature = feature_row$feature,
        total_count = feature_row$total_count,
        samples_present = feature_row$samples_present,
        contrast = info$label,
        reference = info$reference,
        case = info$case,
        samples_used = nrow(current),
        reference_mean_fraction = mean(reference_fraction),
        case_mean_fraction = mean(case_fraction),
        raw_mean_fraction_difference = mean(case_fraction) - mean(reference_fraction),
        log2_odds_ratio = coefficient_table[coefficient, "Estimate"] / log(2),
        adjusted_standard_error_log_odds = coefficient_table[coefficient, "Std. Error"],
        adjusted_p_value = coefficient_table[coefficient, "Pr(>|t|)"],
        stringsAsFactors = FALSE
      )
    }
  }
  comparisons <- if (length(rows)) do.call(rbind, rows) else data.frame()
  if (nrow(comparisons)) {
    comparisons <- add_grouped_fdr(
      comparisons,
      c("feature_type", "receptor", "locus", "contrast")
    )
    comparisons <- comparisons[order(comparisons$contrast, comparisons$FDR_BH), ]
  }
  list(comparisons = comparisons, features = feature_summary)
}

v_usage <- analyse_segment_usage(productive, "v_gene", "V_gene", min_total = 30L, min_samples = 4L)
j_usage <- analyse_segment_usage(productive, "j_gene", "J_gene", min_total = 30L, min_samples = 4L)
vj_usage <- analyse_segment_usage(productive, c("v_gene", "j_gene"), "VJ_pair", min_total = 20L, min_samples = 4L)

write.csv(v_usage$comparisons, file.path(out_dir, "V_gene_clinical_comparisons.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(j_usage$comparisons, file.path(out_dir, "J_gene_clinical_comparisons.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(vj_usage$comparisons, file.path(out_dir, "VJ_pair_clinical_comparisons.csv"), row.names = FALSE, fileEncoding = "UTF-8")
feature_filter_summary <- rbind(v_usage$features, j_usage$features, vj_usage$features)
write.csv(feature_filter_summary, file.path(out_dir, "segment_feature_filter_summary.csv"), row.names = FALSE, fileEncoding = "UTF-8")

linked <- read_tsv_gz(linked_path)
metadata_columns <- intersect(
  c("sample_id", "raw_barcode", "celltype_broad_rechecked", "previous_celltype_main"),
  names(linked)
)
linked_meta <- unique(linked[, metadata_columns, drop = FALSE])
celltype_field <- if ("celltype_broad_rechecked" %in% names(linked_meta)) {
  "celltype_broad_rechecked"
} else {
  "previous_celltype_main"
}
linked_meta$celltype_label <- linked_meta[[celltype_field]]
linked_meta$celltype_label[is.na(linked_meta$celltype_label) | !nzchar(linked_meta$celltype_label)] <- "Unknown"
linked_meta <- linked_meta[, c("sample_id", "raw_barcode", "celltype_label")]
cell_with_type <- merge(
  cell_receptor, linked_meta,
  by = c("sample_id", "raw_barcode"), all.x = TRUE, sort = FALSE
)
cell_with_type$celltype_label[is.na(cell_with_type$celltype_label)] <- "Unknown"

celltype_sample <- aggregate(
  cbind(
    receptor_cells = rep(1L, nrow(cell_with_type)),
    expanded_cells_ge2 = as.integer(cell_with_type$clone_size >= 2),
    expanded_cells_ge5 = as.integer(cell_with_type$clone_size >= 5),
    paired_cells = as.integer(cell_with_type$paired_chain)
  ),
  by = list(
    sample_id = cell_with_type$sample_id,
    receptor = cell_with_type$receptor,
    celltype_label = cell_with_type$celltype_label
  ),
  FUN = sum
)
celltype_sample$expanded_fraction_ge2 <- celltype_sample$expanded_cells_ge2 / celltype_sample$receptor_cells
celltype_sample$expanded_fraction_ge5 <- celltype_sample$expanded_cells_ge5 / celltype_sample$receptor_cells
celltype_sample$paired_fraction <- celltype_sample$paired_cells / celltype_sample$receptor_cells
celltype_sample <- merge(celltype_sample, clinical, by = "sample_id", all.x = TRUE, sort = FALSE)
celltype_sample <- celltype_sample[order(celltype_sample$receptor, celltype_sample$celltype_label, celltype_sample$sample_id), ]
write.csv(celltype_sample, file.path(out_dir, "clone_expansion_by_sample_celltype.csv"), row.names = FALSE, fileEncoding = "UTF-8")

fit_celltype_binomial <- function(data, event_col, outcome_label) {
  groups <- split(data, interaction(data$receptor, data$celltype_label, drop = TRUE, lex.order = TRUE))
  rows <- list()
  index <- 0L
  for (current in groups) {
    current <- current[current$receptor_cells >= 20L, , drop = FALSE]
    if (nrow(current) < 8L) next
    current$events <- current[[event_col]]
    fit <- tryCatch(
      suppressWarnings(glm(
        cbind(events, receptor_cells - events) ~
          batch_wave + severity + allergy_group + plastic_group,
        family = quasibinomial(), data = current
      )),
      error = function(e) NULL
    )
    if (is.null(fit) || fit$rank != length(coef(fit))) next
    coefficient_table <- summary(fit)$coefficients
    current$fraction <- current$events / current$receptor_cells
    for (contrast_name in names(contrast_info)) {
      info <- contrast_info[[contrast_name]]
      coefficient <- info$coefficient
      if (!coefficient %in% rownames(coefficient_table)) next
      group_values <- as.character(current[[info$variable]])
      ref <- current$fraction[group_values == info$reference]
      case <- current$fraction[group_values == info$case]
      index <- index + 1L
      rows[[index]] <- data.frame(
        receptor = unique(current$receptor),
        celltype_label = unique(current$celltype_label),
        outcome = outcome_label,
        contrast = info$label,
        reference = info$reference,
        case = info$case,
        samples_used = nrow(current),
        min_receptor_cells_per_sample = min(current$receptor_cells),
        reference_median_fraction = median(ref),
        case_median_fraction = median(case),
        log2_odds_ratio = coefficient_table[coefficient, "Estimate"] / log(2),
        adjusted_standard_error_log_odds = coefficient_table[coefficient, "Std. Error"],
        adjusted_p_value = coefficient_table[coefficient, "Pr(>|t|)"],
        stringsAsFactors = FALSE
      )
    }
  }
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

celltype_comparisons <- rbind(
  fit_celltype_binomial(celltype_sample, "expanded_cells_ge2", "expanded_fraction_ge2"),
  fit_celltype_binomial(celltype_sample, "expanded_cells_ge5", "expanded_fraction_ge5"),
  fit_celltype_binomial(celltype_sample, "paired_cells", "paired_fraction")
)
if (nrow(celltype_comparisons)) {
  celltype_comparisons <- add_grouped_fdr(celltype_comparisons, c("outcome", "contrast"))
  celltype_comparisons <- celltype_comparisons[order(celltype_comparisons$contrast, celltype_comparisons$FDR_BH), ]
}
write.csv(celltype_comparisons, file.path(out_dir, "celltype_clone_clinical_comparisons.csv"), row.names = FALSE, fileEncoding = "UTF-8")

cdr3_presence <- unique(productive[, c("sample_id", "receptor", "locus", "v_gene", "j_gene", "cdr3_aa")])
cdr3_public <- aggregate(
  sample_id ~ receptor + locus + v_gene + j_gene + cdr3_aa,
  data = cdr3_presence,
  FUN = function(x) length(unique(x))
)
names(cdr3_public)[6] <- "samples_present"
cdr3_counts <- aggregate(
  rep(1L, nrow(productive)),
  by = list(
    receptor = productive$receptor,
    locus = productive$locus,
    v_gene = productive$v_gene,
    j_gene = productive$j_gene,
    cdr3_aa = productive$cdr3_aa
  ),
  FUN = sum
)
names(cdr3_counts)[6] <- "total_contigs"
cdr3_public <- merge(cdr3_public, cdr3_counts, by = c("receptor", "locus", "v_gene", "j_gene", "cdr3_aa"))
cdr3_public <- cdr3_public[cdr3_public$samples_present >= 2L, ]

if (nrow(cdr3_public)) {
  public_keys <- do.call(paste, c(cdr3_public[c("receptor", "locus", "v_gene", "j_gene", "cdr3_aa")], sep = "||"))
  presence_keys <- do.call(paste, c(cdr3_presence[c("receptor", "locus", "v_gene", "j_gene", "cdr3_aa")], sep = "||"))
  cdr3_presence$public_key <- presence_keys
  presence_clinical <- merge(cdr3_presence, clinical, by = "sample_id", all.x = TRUE, sort = FALSE)
  group_count <- function(key, variable, level) {
    sum(presence_clinical$public_key == key & as.character(presence_clinical[[variable]]) == level)
  }
  cdr3_public$轻症_samples <- vapply(public_keys, group_count, integer(1), "severity", "轻症")
  cdr3_public$重症_samples <- vapply(public_keys, group_count, integer(1), "severity", "重症")
  cdr3_public$不过敏_samples <- vapply(public_keys, group_count, integer(1), "allergy_group", "否")
  cdr3_public$过敏_samples <- vapply(public_keys, group_count, integer(1), "allergy_group", "是")
  cdr3_public$非塑形_samples <- vapply(public_keys, group_count, integer(1), "plastic_group", "否")
  cdr3_public$塑形_samples <- vapply(public_keys, group_count, integer(1), "plastic_group", "是")
  cdr3_public <- cdr3_public[order(-cdr3_public$samples_present, -cdr3_public$total_contigs), ]
}
write.csv(cdr3_public, file.path(out_dir, "public_cdr3_vj_sequences.csv"), row.names = FALSE, fileEncoding = "UTF-8")

signature_presence <- unique(cell_receptor[, c("sample_id", "receptor", "clone_signature")])
public_signatures <- aggregate(
  sample_id ~ receptor + clone_signature,
  data = signature_presence,
  FUN = function(x) length(unique(x))
)
names(public_signatures)[3] <- "samples_present"
public_signatures <- public_signatures[public_signatures$samples_present >= 2L, ]
public_signatures <- public_signatures[order(-public_signatures$samples_present, public_signatures$receptor), ]
write.csv(public_signatures, file.path(out_dir, "public_receptor_signatures.csv"), row.names = FALSE, fileEncoding = "UTF-8")

key_metrics <- c(
  "shannon", "inv_simpson", "clonality", "max_clone_fraction",
  "expanded_fraction_ge5", "paired_fraction"
)
metric_plot_data <- do.call(rbind, lapply(key_metrics, function(metric) {
  data.frame(
    sample_id = metrics$sample_id,
    receptor = metrics$receptor,
    metric = metric,
    value = metrics[[metric]],
    severity = as.character(metrics$severity),
    allergy_group = as.character(metrics$allergy_group),
    plastic_group = as.character(metrics$plastic_group),
    stringsAsFactors = FALSE
  )
}))
for (contrast_name in names(contrast_info)) {
  info <- contrast_info[[contrast_name]]
  metric_plot_data$plot_group <- if (contrast_name == "severity") {
    ifelse(metric_plot_data[[info$variable]] == "轻症", "Mild", "Severe")
  } else if (contrast_name == "allergy") {
    ifelse(metric_plot_data[[info$variable]] == "否", "No allergy", "Allergy")
  } else {
    ifelse(metric_plot_data[[info$variable]] == "否", "Non-plastic", "Plastic")
  }
  plot_title <- c(
    severity = "Severe vs mild",
    allergy = "Allergy vs no allergy",
    plastic = "Plastic vs non-plastic"
  )[[contrast_name]]
  plot <- ggplot(metric_plot_data, aes(x = plot_group, y = value, fill = plot_group)) +
    geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.75) +
    geom_jitter(width = 0.12, size = 1.2, alpha = 0.75) +
    facet_grid(metric ~ receptor, scales = "free_y") +
    labs(x = NULL, y = NULL, title = paste("Sample-level repertoire metrics:", plot_title)) +
    theme_bw(base_size = 10) +
    theme(legend.position = "none", strip.text.y = element_text(angle = 0))
  stem <- paste0("repertoire_metrics_", contrast_name)
  ggsave(file.path(out_dir, paste0(stem, ".pdf")), plot, width = 9, height = 13)
  ggsave(file.path(out_dir, paste0(stem, ".png")), plot, width = 9, height = 13, dpi = 240)
}

all_usage <- rbind(v_usage$comparisons, j_usage$comparisons, vj_usage$comparisons)
if (nrow(all_usage)) {
  all_usage$minus_log10_FDR <- -log10(pmax(all_usage$FDR_BH, .Machine$double.xmin))
  all_usage$contrast_plot <- factor(
    all_usage$contrast,
    levels = c("重症_vs_轻症", "过敏_vs_不过敏", "塑形_vs_非塑形"),
    labels = c("Severe vs mild", "Allergy vs no allergy", "Plastic vs non-plastic")
  )
  usage_plot <- ggplot(all_usage, aes(x = log2_odds_ratio, y = minus_log10_FDR, colour = FDR_BH < 0.05)) +
    geom_point(alpha = 0.65, size = 1.2) +
    geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "grey45") +
    facet_grid(feature_type ~ contrast_plot, scales = "free") +
    scale_colour_manual(values = c(`FALSE` = "grey65", `TRUE` = "#B33A3A")) +
    labs(x = "Adjusted log2 odds ratio", y = "-log10(BH FDR)", colour = "FDR < 0.05") +
    theme_bw(base_size = 10)
  ggsave(file.path(out_dir, "segment_usage_clinical_comparisons.pdf"), usage_plot, width = 13, height = 9)
  ggsave(file.path(out_dir, "segment_usage_clinical_comparisons.png"), usage_plot, width = 13, height = 9, dpi = 240)
}

group_counts <- c(
  paste0("轻症 ", sum(clinical$severity == "轻症"), "，重症 ", sum(clinical$severity == "重症")),
  paste0("不过敏 ", sum(clinical$allergy_group == "否"), "，过敏 ", sum(clinical$allergy_group == "是")),
  paste0("非塑形 ", sum(clinical$plastic_group == "否"), "，塑形 ", sum(clinical$plastic_group == "是"))
)

count_significant <- function(data) {
  if (!nrow(data) || !"FDR_BH" %in% names(data)) return(0L)
  sum(data$FDR_BH < 0.05, na.rm = TRUE)
}

report_lines <- c(
  "# 23例 BCR/TCR 临床分组比较",
  "",
  paste0("生成时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## 分析单位与设计",
  "- 所有临床比较都以患儿样本为统计单位，不把单个细胞当作独立患儿。",
  paste0("- 共 23 例；", paste(group_counts, collapse = "；"), "。"),
  paste0("- 校正模型：`~ batch_wave + severity + allergy_group + plastic_group`；设计矩阵秩 ", design_rank, "/", ncol(design), "，满秩。"),
  "- 调整后效应方向均为病例组减参考组：重症-轻症、过敏-不过敏、塑形-非塑形。",
  "- 多重检验使用 Benjamini-Hochberg FDR；同时保留未校正 Wilcoxon P 值用于描述，不以其单独下结论。",
  "",
  "## 已完成内容",
  paste0("- 克隆扩增/多样性：", length(metric_transformations), " 个样本级指标 × BCR/TCR × 3 个临床对比；FDR<0.05 共 ", count_significant(metric_comparisons), " 项。"),
  paste0("- CDR3 长度：按 BCR/TCR 和 IG/TR locus 做样本级均值、中位数、P90；FDR<0.05 共 ", count_significant(cdr3_comparisons), " 项。"),
  paste0("- V 基因使用：过滤后 ", sum(v_usage$features$passed_filter), " 个 receptor-locus-feature；FDR<0.05 共 ", count_significant(v_usage$comparisons), " 项。"),
  paste0("- J 基因使用：过滤后 ", sum(j_usage$features$passed_filter), " 个 receptor-locus-feature；FDR<0.05 共 ", count_significant(j_usage$comparisons), " 项。"),
  paste0("- VJ 配对使用：过滤后 ", sum(vj_usage$features$passed_filter), " 个 receptor-locus-feature；FDR<0.05 共 ", count_significant(vj_usage$comparisons), " 项。"),
  paste0("- 克隆-细胞类型关系：只对每例至少 20 个受体细胞、且至少 8 例可用的组合建模；FDR<0.05 共 ", count_significant(celltype_comparisons), " 项。"),
  paste0("- 跨患儿公共 CDR3-VJ 序列 ", nrow(cdr3_public), " 条；严格整套受体 signature ", nrow(public_signatures), " 条。"),
  "",
  "## 解释限制",
  "- AIRR 数据经审计为混合 VDJ 结果；BCR/TCR 目录多数是重复拷贝，本分析继续按 `locus` 字段拆分 BCR/TCR。",
  "- V/J/VJ 使用采用每例、每个 locus 的计数与总数做 quasibinomial 模型，并经过最低总计数和至少 4 例出现的过滤。",
  "- 公共 CDR3 仅作描述性清单；感染相关公共序列需要独立队列或已知抗原数据库验证。",
  "- 受体出现在非预期 RNA 细胞类型中可能来自双细胞、环境 RNA、细胞类型误差或条码连接误差，不应直接解释为跨谱系生物学。",
  "- 本次 RNA/VDJ 分析不使用旧供应商 WNN 标签；当前矩阵的 RNA+ADT/WNN 已由独立流程重建。"
)
writeLines(report_lines, file.path(out_dir, "BCR_TCR临床分组比较报告.md"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
message("VDJ clinical comparisons completed")
