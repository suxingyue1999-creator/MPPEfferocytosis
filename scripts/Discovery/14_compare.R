#!/usr/bin/env Rscript

lib <- Sys.getenv("BD_R_LIB", "inputs/discovery/bd_rhapsody_10/downstream/10samples_scrna/env/Rlib")
.libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages({
  library(Matrix); library(Seurat); library(edgeR); library(limma)
  library(AnnotationDbi); library(org.Hs.eg.db); library(ggplot2)
})
set.seed(20260726)

root <- Sys.getenv("BD_ANALYSIS_ROOT", "inputs/discovery/bd_rhapsody_10/downstream/10samples_scrna")
result <- file.path(root, "results", "04_compare")
de_dir <- file.path(result, "pseudobulk_de")
enrich_dir <- file.path(result, "enrichment")
plot_dir <- file.path(result, "plots")
dir.create(de_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(enrich_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
object <- readRDS(file.path(root, "results", "03_annotate", "10samples_annotated_seurat_v5.rds"))
DefaultAssay(object) <- "RNA"
object <- tryCatch(JoinLayers(object, assay = "RNA"), error = function(e) object)
counts <- tryCatch(LayerData(object, assay = "RNA", layer = "counts"), error = function(e) NULL)
if (is.null(counts)) counts <- GetAssayData(object, assay = "RNA", slot = "counts")
counts <- as(counts, "dgCMatrix")
meta <- object[[]]
meta$sample_id <- factor(meta$sample_id, levels = c("D01", "D02", "D03", "D04", "D05", "D06", "D07", "D08", "D09", "D10"))
meta$severity <- factor(meta$severity, levels = c("mild", "severe"))
sensitivity_plan <- read.csv(file.path(root, "results", "01_qc", "qc_sensitivity_analysis_plan.csv"), stringsAsFactors = FALSE)
expected_sensitivity_scenarios <- c("sensitivity_exclude_NS4", "sensitivity_exclude_all_qc_outliers")
sensitivity_scenarios <- sensitivity_plan[sensitivity_plan$scenario_id %in% expected_sensitivity_scenarios, , drop = FALSE]
sensitivity_scenarios <- sensitivity_scenarios[match(expected_sensitivity_scenarios, sensitivity_scenarios$scenario_id), , drop = FALSE]
if (nrow(sensitivity_scenarios) != 2L || anyNA(sensitivity_scenarios$scenario_id)) stop("Missing prespecified QC sensitivity scenarios")

sample_severity <- unique(meta[, c("sample_id", "severity"), drop = FALSE])
composition <- as.data.frame(table(sample_id = meta$sample_id, cell_type_broad = meta$cell_type_broad), stringsAsFactors = FALSE)
composition$severity <- sample_severity$severity[match(composition$sample_id, sample_severity$sample_id)]
sample_totals <- aggregate(Freq ~ sample_id, composition, sum)
composition$sample_cells <- sample_totals$Freq[match(composition$sample_id, sample_totals$sample_id)]
composition$proportion <- composition$Freq / composition$sample_cells
write.csv(composition, file.path(result, "sample_cell_type_composition.csv"), row.names = FALSE)

bootstrap_difference <- function(mild, severe, iterations = 5000L) {
  differences <- replicate(iterations, mean(sample(severe, replace = TRUE)) - mean(sample(mild, replace = TRUE)))
  unname(stats::quantile(differences, c(0.025, 0.975), na.rm = TRUE))
}
composition_tests <- lapply(unique(composition$cell_type_broad), function(cell_type) {
  frame <- composition[composition$cell_type_broad == cell_type, , drop = FALSE]
  mild <- frame$proportion[frame$severity == "mild"]
  severe <- frame$proportion[frame$severity == "severe"]
  interval <- bootstrap_difference(mild, severe)
  test <- suppressWarnings(wilcox.test(severe, mild, exact = FALSE, conf.int = FALSE))
  ranks <- rank(c(severe, mild))
  u <- sum(ranks[seq_along(severe)]) - length(severe) * (length(severe) + 1) / 2
  data.frame(
    cell_type_broad = cell_type,
    mild_mean = mean(mild), severe_mean = mean(severe),
    mean_difference_severe_minus_mild = mean(severe) - mean(mild),
    median_difference_severe_minus_mild = median(severe) - median(mild),
    bootstrap_ci_low = interval[1], bootstrap_ci_high = interval[2],
    rank_biserial = 2 * u / (length(severe) * length(mild)) - 1,
    wilcoxon_p_value = test$p.value, stringsAsFactors = FALSE
  )
})
composition_tests <- do.call(rbind, composition_tests)
composition_tests$fdr <- p.adjust(composition_tests$wilcoxon_p_value, method = "BH")
composition_tests <- composition_tests[order(composition_tests$fdr, -abs(composition_tests$mean_difference_severe_minus_mild)), ]
write.csv(composition_tests, file.path(result, "cell_type_composition_comparison.csv"), row.names = FALSE)

p_comp <- ggplot(composition, aes(sample_id, proportion, fill = severity)) +
  geom_col() + facet_wrap(~cell_type_broad, scales = "free_y") +
  scale_fill_manual(values = c(mild = "#2878B5", severe = "#D64B40")) +
  theme_bw(base_size = 10) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(plot_dir, "cell_type_composition_by_patient.pdf"), p_comp, width = 15, height = 10)

cell_count_table <- as.data.frame(table(cell_type_broad = meta$cell_type_broad, sample_id = meta$sample_id), stringsAsFactors = FALSE)
cell_count_table$severity <- sample_severity$severity[match(cell_count_table$sample_id, sample_severity$sample_id)]
cell_count_table <- cell_count_table[cell_count_table$Freq > 0, , drop = FALSE]
write.csv(cell_count_table, file.path(result, "pseudobulk_sample_cell_counts.csv"), row.names = FALSE)

cell_types <- sort(unique(meta$cell_type_broad))
symbol_to_entrez <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = rownames(counts),
  column = "ENTREZID",
  keytype = "SYMBOL",
  multiVals = "first"
)
go_map_all <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(na.omit(symbol_to_entrez)),
  columns = c("GO", "ONTOLOGY"),
  keytype = "ENTREZID"
)
go_map_all <- go_map_all[!is.na(go_map_all$GO) & go_map_all$ONTOLOGY == "BP", , drop = FALSE]
go_ids <- sort(unique(go_map_all$GO))
go_term_map <- setNames(go_ids, go_ids)
cell_type_workers <- suppressWarnings(as.integer(Sys.getenv("BD_CELLTYPE_WORKERS", "10")))
if (!is.finite(cell_type_workers) || cell_type_workers < 1L) cell_type_workers <- 1L
cell_type_workers <- min(cell_type_workers, length(cell_types))

process_cell_type <- function(cell_type) {
  set.seed(20260726 + match(cell_type, cell_types))
  cells <- rownames(meta)[meta$cell_type_broad == cell_type]
  type_meta <- meta[cells, , drop = FALSE]
  sample_cells <- table(type_meta$sample_id)
  eligible_samples <- names(sample_cells)[sample_cells >= 50]
  sample_design <- unique(type_meta[type_meta$sample_id %in% eligible_samples, c("sample_id", "severity"), drop = FALSE])
  sample_design <- sample_design[order(sample_design$sample_id), , drop = FALSE]
  group_counts <- table(sample_design$severity)
  if (length(group_counts) < 2 || any(group_counts[c("mild", "severe")] < 3)) {
    message("Skipping pseudobulk DE for ", cell_type, ": fewer than three eligible patients per group")
    return(NULL)
  }
  selected_cells <- rownames(type_meta)[type_meta$sample_id %in% sample_design$sample_id]
  selected_meta <- type_meta[selected_cells, , drop = FALSE]
  aggregation <- sparse.model.matrix(~ 0 + sample_id, data = selected_meta)
  colnames(aggregation) <- sub("^sample_id", "", colnames(aggregation))
  pseudobulk <- counts[, selected_cells, drop = FALSE] %*% aggregation
  pseudobulk <- as.matrix(pseudobulk[, as.character(sample_design$sample_id), drop = FALSE])
  storage.mode(pseudobulk) <- "integer"
  count_file <- file.path(de_dir, paste0(gsub("[^A-Za-z0-9]+", "_", cell_type), "_pseudobulk_counts.csv.gz"))
  write.csv(data.frame(gene = rownames(pseudobulk), pseudobulk, check.names = FALSE), gzfile(count_file), row.names = FALSE)

  y <- DGEList(counts = pseudobulk, samples = sample_design, group = sample_design$severity)
  keep <- filterByExpr(y, group = sample_design$severity, min.count = 10, min.total.count = 15)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  design <- model.matrix(~ severity, data = sample_design)
  y <- estimateDisp(y, design, robust = TRUE)
  fit <- glmQLFit(y, design, robust = TRUE)
  test <- glmQLFTest(fit, coef = "severitysevere")
  edge <- topTags(test, n = Inf, sort.by = "none")$table
  edge$gene <- rownames(edge)

  voomed <- voom(y, design, plot = FALSE)
  limma_fit <- eBayes(lmFit(voomed, design), robust = TRUE)
  interval <- topTable(limma_fit, coef = "severitysevere", number = Inf, sort.by = "none", confint = 0.95)
  interval$gene <- rownames(interval)
  interval <- interval[match(edge$gene, interval$gene), , drop = FALSE]
  edge$ci95_low <- interval$CI.L
  edge$ci95_high <- interval$CI.R
  edge$cell_type_broad <- cell_type
  edge$contrast <- "severe_vs_mild"
  edge$n_mild_patients <- unname(group_counts[["mild"]])
  edge$n_severe_patients <- unname(group_counts[["severe"]])
  edge$ci_method <- "limma_voom_empirical_bayes"
  edge <- edge[, c("cell_type_broad", "contrast", "gene", "logFC", "ci95_low", "ci95_high", "logCPM", "F", "PValue", "FDR", "n_mild_patients", "n_severe_patients", "ci_method")]
  edge <- edge[order(edge$FDR, -abs(edge$logFC)), , drop = FALSE]
  de_file <- file.path(de_dir, paste0(gsub("[^A-Za-z0-9]+", "_", cell_type), "_edgeR_QLF.csv.gz"))
  write.csv(edge, gzfile(de_file), row.names = FALSE)

  entrez <- symbol_to_entrez[edge$gene]
  go_map <- go_map_all[go_map_all$ENTREZID %in% unique(na.omit(entrez)), , drop = FALSE]
  gene_positions <- split(seq_along(entrez), entrez)
  sets <- lapply(split(go_map$ENTREZID, go_map$GO), function(ids) unique(unlist(gene_positions[unique(ids)], use.names = FALSE)))
  sets <- sets[vapply(sets, length, integer(1)) >= 15 & vapply(sets, length, integer(1)) <= 500]
  signed_statistic <- sign(edge$logFC) * sqrt(pmax(edge$F, 0))
  names(signed_statistic) <- edge$gene
  camera <- cameraPR(signed_statistic, index = sets, use.ranks = TRUE, sort = FALSE)
  camera$go_id <- rownames(camera)
  camera$term <- unname(go_term_map[camera$go_id])
  camera$term[is.na(camera$term) | !nzchar(camera$term)] <- camera$go_id[is.na(camera$term) | !nzchar(camera$term)]
  camera$term_source <- "GO_ID_offline"
  camera$cell_type_broad <- cell_type
  camera$contrast <- "severe_vs_mild"
  camera <- camera[, c("cell_type_broad", "contrast", "go_id", "term", "term_source", "NGenes", "Direction", "PValue", "FDR")]
  camera <- camera[order(camera$FDR, camera$PValue), , drop = FALSE]
  write.csv(camera, gzfile(file.path(enrich_dir, paste0(gsub("[^A-Za-z0-9]+", "_", cell_type), "_GO_BP_cameraPR.csv.gz"))), row.names = FALSE)

  volcano <- edge
  volcano$significant <- volcano$FDR < 0.05 & abs(volcano$logFC) >= 0.5
  p <- ggplot(volcano, aes(logFC, -log10(pmax(FDR, 1e-300)), color = significant)) +
    geom_point(alpha = 0.55, size = 0.8) + scale_color_manual(values = c(`FALSE` = "#777777", `TRUE` = "#C43C39")) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = 2, linewidth = 0.3) +
    geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.3) +
    theme_bw(base_size = 10) + ggtitle(cell_type) + xlab("log2 fold change: severe vs mild") + ylab("-log10 FDR")
  ggsave(file.path(plot_dir, paste0("volcano_", gsub("[^A-Za-z0-9]+", "_", cell_type), ".pdf")), p, width = 7, height = 6)

  run_sensitivity <- function(plan_row) {
    excluded <- strsplit(plan_row$excluded_samples, ";", fixed = TRUE)[[1]]
    sensitivity_design <- sample_design[!as.character(sample_design$sample_id) %in% excluded, , drop = FALSE]
    sensitivity_groups <- table(factor(sensitivity_design$severity, levels = c("mild", "severe")))
    summary_row <- data.frame(
      cell_type_broad = cell_type,
      scenario_id = plan_row$scenario_id,
      excluded_samples = plan_row$excluded_samples,
      status = "SKIPPED",
      detail = "fewer than three eligible patients in at least one severity group",
      genes_tested = 0L,
      fdr_lt_0_05 = 0L,
      n_mild_patients = unname(sensitivity_groups[["mild"]]),
      n_severe_patients = unname(sensitivity_groups[["severe"]]),
      stringsAsFactors = FALSE
    )
    if (any(sensitivity_groups < as.integer(plan_row$min_patients_per_severity))) {
      return(list(de = NULL, summary = summary_row))
    }

    sensitivity_counts <- pseudobulk[, as.character(sensitivity_design$sample_id), drop = FALSE]
    sensitivity_y <- DGEList(counts = sensitivity_counts, samples = sensitivity_design, group = sensitivity_design$severity)
    sensitivity_keep <- filterByExpr(sensitivity_y, group = sensitivity_design$severity, min.count = 10, min.total.count = 15)
    if (sum(sensitivity_keep) < 2L) {
      summary_row$detail <- "fewer than two genes passed expression filtering"
      return(list(de = NULL, summary = summary_row))
    }
    sensitivity_y <- sensitivity_y[sensitivity_keep, , keep.lib.sizes = FALSE]
    sensitivity_y <- calcNormFactors(sensitivity_y)
    sensitivity_model <- model.matrix(~ severity, data = sensitivity_design)
    sensitivity_y <- estimateDisp(sensitivity_y, sensitivity_model, robust = TRUE)
    sensitivity_fit <- glmQLFit(sensitivity_y, sensitivity_model, robust = TRUE)
    sensitivity_test <- glmQLFTest(sensitivity_fit, coef = "severitysevere")
    sensitivity_edge <- topTags(sensitivity_test, n = Inf, sort.by = "none")$table
    sensitivity_edge$gene <- rownames(sensitivity_edge)

    sensitivity_voom <- voom(sensitivity_y, sensitivity_model, plot = FALSE)
    sensitivity_limma <- eBayes(lmFit(sensitivity_voom, sensitivity_model), robust = TRUE)
    sensitivity_interval <- topTable(sensitivity_limma, coef = "severitysevere", number = Inf, sort.by = "none", confint = 0.95)
    sensitivity_interval$gene <- rownames(sensitivity_interval)
    sensitivity_interval <- sensitivity_interval[match(sensitivity_edge$gene, sensitivity_interval$gene), , drop = FALSE]
    sensitivity_edge$ci95_low <- sensitivity_interval$CI.L
    sensitivity_edge$ci95_high <- sensitivity_interval$CI.R
    sensitivity_edge$cell_type_broad <- cell_type
    sensitivity_edge$scenario_id <- plan_row$scenario_id
    sensitivity_edge$excluded_samples <- plan_row$excluded_samples
    sensitivity_edge$contrast <- "severe_vs_mild"
    sensitivity_edge$n_mild_patients <- unname(sensitivity_groups[["mild"]])
    sensitivity_edge$n_severe_patients <- unname(sensitivity_groups[["severe"]])
    sensitivity_edge$ci_method <- "limma_voom_empirical_bayes"
    sensitivity_edge <- sensitivity_edge[, c(
      "scenario_id", "excluded_samples", "cell_type_broad", "contrast", "gene", "logFC",
      "ci95_low", "ci95_high", "logCPM", "F", "PValue", "FDR",
      "n_mild_patients", "n_severe_patients", "ci_method"
    )]
    sensitivity_edge <- sensitivity_edge[order(sensitivity_edge$FDR, -abs(sensitivity_edge$logFC)), , drop = FALSE]
    sensitivity_file <- file.path(
      de_dir,
      paste0(gsub("[^A-Za-z0-9]+", "_", cell_type), "_", plan_row$scenario_id, "_edgeR_QLF.csv.gz")
    )
    write.csv(sensitivity_edge, gzfile(sensitivity_file), row.names = FALSE)
    summary_row$status <- "COMPLETED"
    summary_row$detail <- "edgeR QL and limma-voom interval completed"
    summary_row$genes_tested <- nrow(sensitivity_edge)
    summary_row$fdr_lt_0_05 <- sum(sensitivity_edge$FDR < 0.05)
    list(de = sensitivity_edge, summary = summary_row)
  }

  sensitivity_runs <- lapply(seq_len(nrow(sensitivity_scenarios)), function(index) {
    run_sensitivity(sensitivity_scenarios[index, , drop = FALSE])
  })
  sensitivity_de <- Filter(Negate(is.null), lapply(sensitivity_runs, `[[`, "de"))
  sensitivity_de <- if (length(sensitivity_de)) do.call(rbind, sensitivity_de) else NULL
  sensitivity_summary <- do.call(rbind, lapply(sensitivity_runs, `[[`, "summary"))
  summary <- data.frame(
    cell_type_broad = cell_type, genes_tested = nrow(edge),
    fdr_lt_0_05 = sum(edge$FDR < 0.05), fdr_lt_0_05_abs_logfc_ge_0_5 = sum(edge$FDR < 0.05 & abs(edge$logFC) >= 0.5),
    n_mild_patients = unname(group_counts[["mild"]]), n_severe_patients = unname(group_counts[["severe"]])
  )
  list(
    cell_type = cell_type,
    de = edge,
    enrichment = camera,
    summary = summary,
    sensitivity_de = sensitivity_de,
    sensitivity_summary = sensitivity_summary
  )
}

message("Pseudobulk cell-type workers: ", cell_type_workers)
process_cell_type_safe <- function(cell_type) {
  tryCatch(
    process_cell_type(cell_type),
    error = function(error) structure(
      list(cell_type = cell_type, error = conditionMessage(error)),
      class = "pseudobulk_worker_error"
    )
  )
}
cell_type_results <- parallel::mclapply(
  cell_types,
  process_cell_type_safe,
  mc.cores = cell_type_workers,
  mc.preschedule = FALSE,
  mc.set.seed = TRUE
)
failed_cell_types <- vapply(cell_type_results, inherits, logical(1), what = "pseudobulk_worker_error")
if (any(failed_cell_types)) {
  failure_report <- data.frame(
    cell_type_broad = vapply(cell_type_results[failed_cell_types], `[[`, character(1), "cell_type"),
    error = vapply(cell_type_results[failed_cell_types], `[[`, character(1), "error"),
    stringsAsFactors = FALSE
  )
  write.csv(failure_report, file.path(result, "pseudobulk_worker_failures.csv"), row.names = FALSE)
  stop("Pseudobulk worker failure: ", paste(failure_report$cell_type_broad, collapse = ","))
}
cell_type_results <- Filter(Negate(is.null), cell_type_results)
all_de <- setNames(lapply(cell_type_results, `[[`, "de"), vapply(cell_type_results, `[[`, character(1), "cell_type"))
all_enrichment <- setNames(lapply(cell_type_results, `[[`, "enrichment"), vapply(cell_type_results, `[[`, character(1), "cell_type"))
de_summary <- setNames(lapply(cell_type_results, `[[`, "summary"), vapply(cell_type_results, `[[`, character(1), "cell_type"))
sensitivity_de <- Filter(Negate(is.null), lapply(cell_type_results, `[[`, "sensitivity_de"))
sensitivity_summary <- lapply(cell_type_results, `[[`, "sensitivity_summary")

if (!length(all_de)) stop("No cell type met pseudobulk patient-replicate guardrails")
combined_de <- do.call(rbind, all_de)
combined_enrichment <- do.call(rbind, all_enrichment)
write.csv(combined_de, gzfile(file.path(result, "all_cell_types_pseudobulk_de.csv.gz")), row.names = FALSE)
write.csv(combined_enrichment, gzfile(file.path(result, "all_cell_types_GO_BP_cameraPR.csv.gz")), row.names = FALSE)
write.csv(do.call(rbind, de_summary), file.path(result, "pseudobulk_de_summary.csv"), row.names = FALSE)
combined_sensitivity_de <- if (length(sensitivity_de)) do.call(rbind, sensitivity_de) else data.frame()
combined_sensitivity_summary <- do.call(rbind, sensitivity_summary)
write.csv(combined_sensitivity_de, gzfile(file.path(result, "all_cell_types_pseudobulk_sensitivity_de.csv.gz")), row.names = FALSE)
write.csv(combined_sensitivity_summary, file.path(result, "pseudobulk_sensitivity_summary.csv"), row.names = FALSE)

sensitivity_concordance <- do.call(rbind, lapply(split(combined_sensitivity_de, list(combined_sensitivity_de$cell_type_broad, combined_sensitivity_de$scenario_id), drop = TRUE), function(frame) {
  primary <- combined_de[combined_de$cell_type_broad == frame$cell_type_broad[[1]], c("gene", "logFC"), drop = FALSE]
  matched <- merge(primary, frame[, c("gene", "logFC"), drop = FALSE], by = "gene", suffixes = c("_primary", "_sensitivity"))
  data.frame(
    cell_type_broad = frame$cell_type_broad[[1]],
    scenario_id = frame$scenario_id[[1]],
    genes_compared = nrow(matched),
    logFC_spearman = suppressWarnings(cor(matched$logFC_primary, matched$logFC_sensitivity, method = "spearman")),
    direction_agreement = mean(sign(matched$logFC_primary) == sign(matched$logFC_sensitivity)),
    stringsAsFactors = FALSE
  )
}))
write.csv(sensitivity_concordance, file.path(result, "pseudobulk_sensitivity_concordance.csv"), row.names = FALSE)

top_pathways <- combined_enrichment[combined_enrichment$FDR < 0.10, , drop = FALSE]
if (nrow(top_pathways)) {
  top_pathways <- do.call(rbind, lapply(split(top_pathways, top_pathways$cell_type_broad), head, 10))
  p_path <- ggplot(top_pathways, aes(-log10(pmax(FDR, 1e-300)), reorder(term, -log10(pmax(FDR, 1e-300))), color = Direction, size = NGenes)) +
    geom_point() + facet_wrap(~cell_type_broad, scales = "free_y") + theme_bw(base_size = 9) + ylab(NULL) + xlab("-log10 FDR")
  ggsave(file.path(plot_dir, "GO_BP_cameraPR_top_pathways.pdf"), p_path, width = 14, height = 10)
} else {
  writeLines("No GO BP pathway reached FDR < 0.10.", file.path(plot_dir, "GO_BP_no_fdr_lt_0.10.txt"))
}

writeLines(c(
  "# Stage 04: mild versus severe comparison",
  "",
  "- The patient sample is the statistical replicate; cells are never treated as independent replicates.",
  paste0("- Parallelism: ", cell_type_workers, " independent broad-cell-type pseudobulk workers."),
  "- Composition effects are severe-minus-mild patient-level proportions with bootstrap 95% confidence intervals, Wilcoxon P values, and BH FDR.",
  "- Pseudobulk matrices sum unmodified RNA counts within each patient and broad cell type.",
  "- A cell type is tested only with at least 50 cells in at least three patients per severity group.",
  "- edgeR quasi-likelihood tests provide primary P values/FDR; limma-voom empirical-Bayes intervals provide 95% log2-fold-change confidence intervals.",
  "- GO Biological Process competitive enrichment uses limma cameraPR and local org.Hs.eg.db membership; the fixed offline environment reports standardized GO identifiers in the term column (term_source=GO_ID_offline).",
  "- Prespecified sensitivity analyses rerun patient-level pseudobulk DE after excluding D04 alone and after excluding D03, D04, and D07 together; each scenario still requires at least three patients per severity group.",
  "- Positive logFC and pathway direction Up indicate higher expression in severe disease.",
  "- With five patients per group, effect sizes and intervals should be emphasized over significance alone."
), file.path(result, "comparison_report.md"))
writeLines(capture.output(sessionInfo()), file.path(result, "r_sessionInfo.txt"))
