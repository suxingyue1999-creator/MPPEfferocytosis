#!/usr/bin/env Rscript

lib <- Sys.getenv("BD_R_LIB", "inputs/discovery/bd_rhapsody_10/downstream/10samples_scrna/env/host-Rlib")
.libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages({
  library(limma)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(ggplot2)
})

options(stringsAsFactors = FALSE, warn = 1)
set.seed(20260727)

cross_root <- Sys.getenv(
  "BD_CROSS_ROOT",
  "inputs/discovery/bd_rhapsody_10/downstream/10_23_cross_cohort_v1"
)
stage02 <- file.path(cross_root, "results", "02_harmonized_validation")
result <- file.path(cross_root, "results", "03_meta_pathway")
figure_root <- file.path(result, "figures")
dir.create(figure_root, recursive = TRUE, showWarnings = FALSE)

read_csv_path <- function(path) {
  if (!file.exists(path)) stop("Missing required input: ", path)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8-BOM")
}

read_gz_csv_path <- function(path) {
  if (!file.exists(path)) stop("Missing required input: ", path)
  read.csv(gzfile(path), stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8-BOM")
}

write_csv <- function(frame, name) {
  write.csv(frame, file.path(result, name), row.names = FALSE, fileEncoding = "UTF-8")
}

primary <- read_gz_csv_path(file.path(stage02, "harmonized_primary_de.csv.gz"))
adjusted <- read_gz_csv_path(file.path(stage02, "validation_adjusted_de.csv.gz"))
composition <- read_csv_path(file.path(stage02, "harmonized_composition_comparison.csv"))
stage02_inventory <- read_csv_path(file.path(stage02, "stage_02_input_inventory.csv"))
lineages <- c("B_lineage", "Macrophage_Monocyte", "Neutrophil", "T_lineage")

if (!all(primary$effect_se > 0 & is.finite(primary$effect_se))) stop("Primary effects contain invalid standard errors")
if (!setequal(unique(primary$harmonized_lineage), lineages)) stop("Primary effect lineages are incomplete")

effect_frame <- function(frame, cohort, prefix) {
  selected <- frame[frame$cohort == cohort, c(
    "harmonized_lineage", "gene", "edgeR_logFC", "effect_log2fc", "effect_se",
    "PValue", "FDR", "F", "n_mild_patients", "n_severe_patients"
  ), drop = FALSE]
  names(selected)[-c(1, 2)] <- paste0(prefix, "_", names(selected)[-c(1, 2)])
  selected
}

discovery <- effect_frame(primary, "discovery_10", "discovery")
validation <- effect_frame(primary, "validation_23", "validation")
gene_meta <- merge(discovery, validation, by = c("harmonized_lineage", "gene"), all = FALSE, sort = FALSE)
if (anyDuplicated(gene_meta[, c("harmonized_lineage", "gene")])) stop("Gene meta input contains duplicate keys")

meta_two_vectorized <- function(effect1, se1, effect2, se2) {
  weight1 <- 1 / se1^2
  weight2 <- 1 / se2^2
  weight_sum <- weight1 + weight2
  fixed_effect <- (weight1 * effect1 + weight2 * effect2) / weight_sum
  fixed_se <- sqrt(1 / weight_sum)
  q <- weight1 * (effect1 - fixed_effect)^2 + weight2 * (effect2 - fixed_effect)^2
  c_value <- weight_sum - (weight1^2 + weight2^2) / weight_sum
  tau2 <- pmax(0, (q - 1) / c_value)
  random_weight1 <- 1 / (se1^2 + tau2)
  random_weight2 <- 1 / (se2^2 + tau2)
  random_weight_sum <- random_weight1 + random_weight2
  random_effect <- (random_weight1 * effect1 + random_weight2 * effect2) / random_weight_sum
  random_se <- sqrt(1 / random_weight_sum)
  data.frame(
    fixed_effect = fixed_effect,
    fixed_se = fixed_se,
    fixed_z = fixed_effect / fixed_se,
    fixed_p_value = 2 * pnorm(-abs(fixed_effect / fixed_se)),
    q = q,
    q_p_value = pchisq(q, df = 1, lower.tail = FALSE),
    i2_percent = ifelse(q > 0, pmax(0, (q - 1) / q) * 100, 0),
    tau2 = tau2,
    random_effect = random_effect,
    random_se = random_se,
    random_z = random_effect / random_se,
    random_p_value = 2 * pnorm(-abs(random_effect / random_se)),
    stringsAsFactors = FALSE
  )
}

meta_values <- meta_two_vectorized(
  gene_meta$discovery_effect_log2fc,
  gene_meta$discovery_effect_se,
  gene_meta$validation_effect_log2fc,
  gene_meta$validation_effect_se
)
gene_meta <- cbind(gene_meta, meta_values)
gene_meta$fixed_ci95_low <- gene_meta$fixed_effect - 1.96 * gene_meta$fixed_se
gene_meta$fixed_ci95_high <- gene_meta$fixed_effect + 1.96 * gene_meta$fixed_se
gene_meta$random_ci95_low <- gene_meta$random_effect - 1.96 * gene_meta$random_se
gene_meta$random_ci95_high <- gene_meta$random_effect + 1.96 * gene_meta$random_se
gene_meta$fixed_fdr <- ave(gene_meta$fixed_p_value, gene_meta$harmonized_lineage, FUN = function(x) p.adjust(x, "BH"))
gene_meta$random_fdr <- ave(gene_meta$random_p_value, gene_meta$harmonized_lineage, FUN = function(x) p.adjust(x, "BH"))
gene_meta$direction_concordant <- sign(gene_meta$discovery_effect_log2fc) == sign(gene_meta$validation_effect_log2fc)
gene_meta$gene_replicated <- with(
  gene_meta,
  direction_concordant & discovery_FDR < 0.05 & validation_PValue < 0.05 & fixed_fdr < 0.05 & i2_percent < 75
)
gene_meta$replication_class <- "discordant_or_unsupported"
gene_meta$replication_class[gene_meta$direction_concordant] <- "directionally_consistent_not_significant"
gene_meta$replication_class[
  gene_meta$direction_concordant & gene_meta$discovery_FDR < 0.05 & gene_meta$validation_PValue >= 0.05
] <- "discovery_only"
gene_meta$replication_class[
  gene_meta$direction_concordant & gene_meta$discovery_PValue < 0.05 & gene_meta$validation_PValue < 0.05 & gene_meta$fixed_fdr < 0.05
] <- "meta_supported_concordant"
gene_meta$replication_class[gene_meta$gene_replicated] <- "replicated_gene"
gene_meta <- gene_meta[order(gene_meta$harmonized_lineage, gene_meta$fixed_fdr, -abs(gene_meta$fixed_effect)), , drop = FALSE]
write.csv(gene_meta, gzfile(file.path(result, "two_cohort_gene_meta_analysis.csv.gz")), row.names = FALSE)

gene_concordance <- do.call(rbind, lapply(split(gene_meta, gene_meta$harmonized_lineage), function(frame) {
  data.frame(
    harmonized_lineage = frame$harmonized_lineage[[1]],
    common_genes = nrow(frame),
    effect_spearman = cor(frame$discovery_effect_log2fc, frame$validation_effect_log2fc, method = "spearman"),
    effect_pearson = cor(frame$discovery_effect_log2fc, frame$validation_effect_log2fc, method = "pearson"),
    direction_agreement = mean(frame$direction_concordant),
    discovery_fdr_hits = sum(frame$discovery_FDR < 0.05),
    validation_fdr_hits = sum(frame$validation_FDR < 0.05),
    fixed_meta_fdr_hits = sum(frame$fixed_fdr < 0.05),
    replicated_genes = sum(frame$gene_replicated),
    stringsAsFactors = FALSE
  )
}))
write_csv(gene_concordance, "gene_effect_concordance_summary.csv")

validation_primary <- effect_frame(primary, "validation_23", "primary")
adjusted_selected <- adjusted[, c("harmonized_lineage", "gene", "effect_log2fc", "effect_se", "PValue", "FDR"), drop = FALSE]
names(adjusted_selected)[-c(1, 2)] <- paste0("adjusted_", names(adjusted_selected)[-c(1, 2)])
adjustment_comparison <- merge(validation_primary, adjusted_selected, by = c("harmonized_lineage", "gene"), all = FALSE)
adjustment_summary <- do.call(rbind, lapply(split(adjustment_comparison, adjustment_comparison$harmonized_lineage), function(frame) {
  data.frame(
    harmonized_lineage = frame$harmonized_lineage[[1]],
    common_genes = nrow(frame),
    effect_spearman = cor(frame$primary_effect_log2fc, frame$adjusted_effect_log2fc, method = "spearman"),
    direction_agreement = mean(sign(frame$primary_effect_log2fc) == sign(frame$adjusted_effect_log2fc)),
    adjusted_nominal_hits = sum(frame$adjusted_PValue < 0.05),
    adjusted_fdr_hits = sum(frame$adjusted_FDR < 0.05),
    stringsAsFactors = FALSE
  )
}))
write_csv(adjustment_summary, "validation_adjustment_concordance.csv")

validation_root <- sub(
  "/results/.*$",
  "",
  stage02_inventory$path[grepl("pseudobulk_raw_counts[.]rds$", stage02_inventory$path)]
)
term_dictionary <- data.frame(go_id = character(), term = character(), stringsAsFactors = FALSE)
if (length(validation_root) == 1L && nzchar(validation_root)) {
  local_go_files <- list.files(
    file.path(validation_root, "results", "08_lineage_gsea"),
    pattern = "GO_BP_GSEA[.]csv$",
    recursive = TRUE,
    full.names = TRUE
  )
  local_terms <- lapply(local_go_files, function(path) {
    frame <- tryCatch(read_csv_path(path), error = function(e) NULL)
    if (is.null(frame) || !all(c("ID", "Description") %in% names(frame))) return(NULL)
    frame[, c("ID", "Description"), drop = FALSE]
  })
  local_terms <- Filter(Negate(is.null), local_terms)
  if (length(local_terms)) {
    term_dictionary <- unique(do.call(rbind, local_terms))
    names(term_dictionary) <- c("go_id", "term")
    term_dictionary <- term_dictionary[!duplicated(term_dictionary$go_id), , drop = FALSE]
  }
}

make_go_sets <- function(genes) {
  symbol_to_entrez <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = genes,
    column = "ENTREZID",
    keytype = "SYMBOL",
    multiVals = "first"
  )
  go_map <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = unique(na.omit(symbol_to_entrez)),
    columns = c("GO", "ONTOLOGY"),
    keytype = "ENTREZID"
  )
  go_map <- go_map[!is.na(go_map$GO) & go_map$ONTOLOGY == "BP", , drop = FALSE]
  gene_positions <- split(seq_along(symbol_to_entrez), symbol_to_entrez)
  sets <- lapply(split(go_map$ENTREZID, go_map$GO), function(ids) {
    unique(unlist(gene_positions[unique(ids)], use.names = FALSE))
  })
  sizes <- vapply(sets, length, integer(1))
  sets[sizes >= 15L & sizes <= 500L]
}

pathway_results <- list()
pathway_meta_results <- list()
pathway_index <- 0L
for (lineage in lineages) {
  discovery_frame <- primary[primary$cohort == "discovery_10" & primary$harmonized_lineage == lineage, , drop = FALSE]
  validation_frame <- primary[primary$cohort == "validation_23" & primary$harmonized_lineage == lineage, , drop = FALSE]
  common_genes <- sort(intersect(discovery_frame$gene, validation_frame$gene))
  discovery_frame <- discovery_frame[match(common_genes, discovery_frame$gene), , drop = FALSE]
  validation_frame <- validation_frame[match(common_genes, validation_frame$gene), , drop = FALSE]
  go_sets <- make_go_sets(common_genes)
  if (length(go_sets) < 100L) stop(lineage, " has too few offline GO-BP sets")

  cohort_camera <- lapply(list(discovery_10 = discovery_frame, validation_23 = validation_frame), function(frame) {
    signed_statistic <- sign(frame$effect_log2fc) * sqrt(pmax(frame$F, 0))
    names(signed_statistic) <- frame$gene
    camera <- cameraPR(signed_statistic, index = go_sets, use.ranks = TRUE, sort = FALSE)
    camera$go_id <- rownames(camera)
    camera$term <- rep(NA_character_, nrow(camera))
    if (nrow(term_dictionary)) {
      camera$term <- term_dictionary$term[match(camera$go_id, term_dictionary$go_id)]
    }
    camera$term[is.na(camera$term) | !nzchar(camera$term)] <- camera$go_id[is.na(camera$term) | !nzchar(camera$term)]
    camera$harmonized_lineage <- lineage
    camera$cohort <- frame$cohort[[1]]
    camera$common_gene_universe <- length(common_genes)
    camera$term_source <- "GO_BP_org.Hs.eg.db_offline_common_universe"
    camera[, c(
      "cohort", "harmonized_lineage", "go_id", "term", "term_source",
      "common_gene_universe", "NGenes", "Direction", "PValue", "FDR"
    )]
  })
  pathway_index <- pathway_index + 1L
  pathway_results[[pathway_index]] <- cohort_camera$discovery_10
  pathway_index <- pathway_index + 1L
  pathway_results[[pathway_index]] <- cohort_camera$validation_23

  d <- cohort_camera$discovery_10
  v <- cohort_camera$validation_23
  v <- v[match(d$go_id, v$go_id), , drop = FALSE]
  if (anyNA(v$go_id)) stop("Shared GO set alignment failed for ", lineage)
  signed_z_discovery <- ifelse(d$Direction == "Up", 1, -1) * qnorm(pmax(d$PValue, 1e-300) / 2, lower.tail = FALSE)
  signed_z_validation <- ifelse(v$Direction == "Up", 1, -1) * qnorm(pmax(v$PValue, 1e-300) / 2, lower.tail = FALSE)
  n_discovery <- unique(discovery_frame$n_mild_patients + discovery_frame$n_severe_patients)
  n_validation <- unique(validation_frame$n_mild_patients + validation_frame$n_severe_patients)
  combined_z <- (sqrt(n_discovery) * signed_z_discovery + sqrt(n_validation) * signed_z_validation) /
    sqrt(n_discovery + n_validation)
  pathway_meta <- data.frame(
    harmonized_lineage = lineage,
    go_id = d$go_id,
    term = d$term,
    NGenes = d$NGenes,
    common_gene_universe = d$common_gene_universe,
    discovery_direction = d$Direction,
    discovery_p_value = d$PValue,
    discovery_fdr = d$FDR,
    validation_direction = v$Direction,
    validation_p_value = v$PValue,
    validation_fdr = v$FDR,
    direction_concordant = d$Direction == v$Direction,
    combined_z = combined_z,
    combined_p_value = 2 * pnorm(-abs(combined_z)),
    stringsAsFactors = FALSE
  )
  pathway_meta$combined_fdr <- p.adjust(pathway_meta$combined_p_value, "BH")
  pathway_meta$strictly_replicated_pathway <- with(
    pathway_meta,
    direction_concordant & discovery_fdr < 0.05 & validation_fdr < 0.05
  )
  pathway_meta$supported_pathway <- with(
    pathway_meta,
    direction_concordant & discovery_fdr < 0.05 & validation_p_value < 0.05 & combined_fdr < 0.05
  )
  pathway_meta_results[[lineage]] <- pathway_meta
}

pathway_all <- do.call(rbind, pathway_results)
pathway_meta <- do.call(rbind, pathway_meta_results)
rownames(pathway_meta) <- NULL
write.csv(pathway_all, gzfile(file.path(result, "shared_GO_BP_cameraPR.csv.gz")), row.names = FALSE)
write.csv(pathway_meta, gzfile(file.path(result, "two_cohort_pathway_integration.csv.gz")), row.names = FALSE)
write_csv(term_dictionary, "local_GO_term_dictionary.csv")

replication_summary <- merge(gene_concordance, do.call(rbind, lapply(split(pathway_meta, pathway_meta$harmonized_lineage), function(frame) {
  data.frame(
    harmonized_lineage = frame$harmonized_lineage[[1]],
    shared_pathways = nrow(frame),
    strictly_replicated_pathways = sum(frame$strictly_replicated_pathway),
    supported_pathways = sum(frame$supported_pathway),
    stringsAsFactors = FALSE
  )
})), by = "harmonized_lineage", all = TRUE)
replication_summary$mechanism_eligible <- with(
  replication_summary,
  replicated_genes >= 3L | strictly_replicated_pathways >= 1L
)
replication_summary$replication_status <- "not_replicated"
replication_summary$replication_status[
  !replication_summary$mechanism_eligible &
    (replication_summary$replicated_genes > 0L | replication_summary$supported_pathways > 0L)
] <- "replicated_limited"
replication_summary$replication_status[replication_summary$mechanism_eligible] <- "replicated_strong"
replication_summary$decision_rule <- paste(
  "mechanism eligible if replicated_genes >= 3 or strictly_replicated_pathways >= 1;",
  "replicated gene requires discovery FDR < 0.05, validation P < 0.05, concordant direction,",
  "fixed-effect meta FDR < 0.05, and I2 < 75%; strict pathway requires concordant direction and FDR < 0.05 in both cohorts"
)
write_csv(replication_summary, "lineage_replication_decisions.csv")

composition_d <- composition[composition$cohort == "discovery_10", ]
composition_v <- composition[composition$cohort == "validation_23", ]
composition_meta <- merge(
  composition_d[, c("harmonized_lineage", "mean_difference", "effect_se", "wilcoxon_p_value", "fdr_within_cohort")],
  composition_v[, c("harmonized_lineage", "mean_difference", "effect_se", "wilcoxon_p_value", "fdr_within_cohort")],
  by = "harmonized_lineage",
  suffixes = c("_discovery", "_validation")
)
composition_values <- meta_two_vectorized(
  composition_meta$mean_difference_discovery,
  composition_meta$effect_se_discovery,
  composition_meta$mean_difference_validation,
  composition_meta$effect_se_validation
)
composition_meta <- cbind(composition_meta, composition_values)
composition_meta$fixed_fdr <- p.adjust(composition_meta$fixed_p_value, "BH")
composition_meta$random_fdr <- p.adjust(composition_meta$random_p_value, "BH")
composition_meta$direction_concordant <- sign(composition_meta$mean_difference_discovery) == sign(composition_meta$mean_difference_validation)
write_csv(composition_meta, "two_cohort_composition_meta_analysis.csv")

scatter_plot <- ggplot(gene_meta, aes(discovery_effect_log2fc, validation_effect_log2fc)) +
  geom_hline(yintercept = 0, color = "#B7B7B7", linewidth = 0.25) +
  geom_vline(xintercept = 0, color = "#B7B7B7", linewidth = 0.25) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "#333333", linewidth = 0.35) +
  geom_point(aes(color = gene_replicated), alpha = 0.35, size = 0.55) +
  facet_wrap(~ harmonized_lineage, scales = "free") +
  scale_color_manual(values = c(`FALSE` = "#7A7A7A", `TRUE` = "#C43C39")) +
  theme_bw(base_size = 10) +
  labs(x = "Discovery log2 effect", y = "Validation log2 effect", color = "Replicated gene")
ggsave(file.path(figure_root, "gene_effect_concordance.pdf"), scatter_plot, width = 11, height = 8, useDingbats = FALSE)
ggsave(file.path(figure_root, "gene_effect_concordance.png"), scatter_plot, width = 11, height = 8, dpi = 180)

top_gene_rows <- do.call(rbind, lapply(split(gene_meta, gene_meta$harmonized_lineage), function(frame) head(frame[order(frame$fixed_fdr), ], 6L)))
forest_frame <- rbind(
  data.frame(harmonized_lineage = top_gene_rows$harmonized_lineage, gene = top_gene_rows$gene, estimate = "Discovery", effect = top_gene_rows$discovery_effect_log2fc, low = top_gene_rows$discovery_effect_log2fc - 1.96 * top_gene_rows$discovery_effect_se, high = top_gene_rows$discovery_effect_log2fc + 1.96 * top_gene_rows$discovery_effect_se),
  data.frame(harmonized_lineage = top_gene_rows$harmonized_lineage, gene = top_gene_rows$gene, estimate = "Validation", effect = top_gene_rows$validation_effect_log2fc, low = top_gene_rows$validation_effect_log2fc - 1.96 * top_gene_rows$validation_effect_se, high = top_gene_rows$validation_effect_log2fc + 1.96 * top_gene_rows$validation_effect_se),
  data.frame(harmonized_lineage = top_gene_rows$harmonized_lineage, gene = top_gene_rows$gene, estimate = "Fixed meta", effect = top_gene_rows$fixed_effect, low = top_gene_rows$fixed_ci95_low, high = top_gene_rows$fixed_ci95_high)
)
forest_frame$gene <- factor(forest_frame$gene, levels = rev(unique(top_gene_rows$gene)))
forest_plot <- ggplot(forest_frame, aes(gene, effect, color = estimate)) +
  geom_hline(yintercept = 0, color = "#777777", linewidth = 0.3) +
  geom_errorbar(aes(ymin = low, ymax = high), position = position_dodge(width = 0.55), width = 0.2) +
  geom_point(position = position_dodge(width = 0.55), size = 1.4) +
  coord_flip() +
  facet_wrap(~ harmonized_lineage, scales = "free_y") +
  scale_color_manual(values = c("Discovery" = "#2878B5", "Validation" = "#D64B40", "Fixed meta" = "#222222")) +
  theme_bw(base_size = 10) +
  labs(x = NULL, y = "log2 effect with 95% CI", color = NULL)
ggsave(file.path(figure_root, "top_gene_meta_forest.pdf"), forest_plot, width = 12, height = 9, useDingbats = FALSE)
ggsave(file.path(figure_root, "top_gene_meta_forest.png"), forest_plot, width = 12, height = 9, dpi = 180)

top_pathways <- do.call(rbind, lapply(split(pathway_meta, pathway_meta$harmonized_lineage), function(frame) head(frame[order(frame$combined_fdr), ], 8L)))
top_pathways$display_term <- ifelse(nchar(top_pathways$term) > 60, paste0(substr(top_pathways$term, 1, 57), "..."), top_pathways$term)
top_pathways$signed_strength <- sign(top_pathways$combined_z) * -log10(pmax(top_pathways$combined_fdr, 1e-300))
pathway_plot <- ggplot(top_pathways, aes(reorder(display_term, signed_strength), signed_strength, fill = strictly_replicated_pathway)) +
  geom_col(width = 0.75) +
  coord_flip() +
  facet_wrap(~ harmonized_lineage, scales = "free_y") +
  scale_fill_manual(values = c(`FALSE` = "#6E7781", `TRUE` = "#C43C39")) +
  theme_bw(base_size = 9) +
  labs(x = NULL, y = "Signed -log10 combined FDR", fill = "Strict replication")
ggsave(file.path(figure_root, "shared_pathway_integration.pdf"), pathway_plot, width = 13, height = 10, useDingbats = FALSE)
ggsave(file.path(figure_root, "shared_pathway_integration.png"), pathway_plot, width = 13, height = 10, dpi = 180)

composition_forest <- rbind(
  data.frame(harmonized_lineage = composition_meta$harmonized_lineage, estimate = "Discovery", effect = composition_meta$mean_difference_discovery, se = composition_meta$effect_se_discovery),
  data.frame(harmonized_lineage = composition_meta$harmonized_lineage, estimate = "Validation", effect = composition_meta$mean_difference_validation, se = composition_meta$effect_se_validation),
  data.frame(harmonized_lineage = composition_meta$harmonized_lineage, estimate = "Fixed meta", effect = composition_meta$fixed_effect, se = composition_meta$fixed_se)
)
composition_forest$low <- composition_forest$effect - 1.96 * composition_forest$se
composition_forest$high <- composition_forest$effect + 1.96 * composition_forest$se
composition_plot <- ggplot(composition_forest, aes(harmonized_lineage, effect, color = estimate)) +
  geom_hline(yintercept = 0, color = "#777777", linewidth = 0.3) +
  geom_errorbar(aes(ymin = low, ymax = high), position = position_dodge(width = 0.55), width = 0.18) +
  geom_point(position = position_dodge(width = 0.55), size = 1.7) +
  coord_flip() +
  scale_color_manual(values = c("Discovery" = "#2878B5", "Validation" = "#D64B40", "Fixed meta" = "#222222")) +
  theme_bw(base_size = 10) +
  labs(x = NULL, y = "Severe minus mild cell proportion", color = NULL)
ggsave(file.path(figure_root, "composition_meta_forest.pdf"), composition_plot, width = 9, height = 5.5, useDingbats = FALSE)
ggsave(file.path(figure_root, "composition_meta_forest.png"), composition_plot, width = 9, height = 5.5, dpi = 180)

eligible <- replication_summary$harmonized_lineage[replication_summary$mechanism_eligible]
report <- c(
  "# Cross-cohort effect and pathway integration",
  "",
  paste0("- Generated: ", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
  "- Cohorts are independent; two-cohort fixed and DerSimonian-Laird random effects are reported.",
  "- Gene effects use limma-voom coefficients and moderated standard errors from the harmonized patient-level models.",
  "- Differential-expression evidence remains the edgeR robust QL P value and FDR.",
  "- GO-BP testing uses cameraPR with ranks on the exact shared tested-gene universe for each lineage.",
  "- GO membership comes from the local org.Hs.eg.db; local prior GSEA outputs are used only as an optional GO-ID description dictionary.",
  "- Gene replication: discovery FDR < 0.05, validation P < 0.05, concordant direction, fixed meta FDR < 0.05, and I2 < 75%.",
  "- Strict pathway replication: concordant direction and cameraPR FDR < 0.05 in both cohorts.",
  "- Mechanism eligibility: at least three replicated genes or at least one strictly replicated pathway.",
  paste0("- Mechanism-eligible lineages: ", if (length(eligible)) paste(eligible, collapse = ", ") else "none"),
  "",
  "A pooled meta-analysis result alone is not labeled independent replication when validation support or direction concordance is absent."
)
writeLines(report, file.path(result, "meta_pathway_report.md"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(result, "sessionInfo.txt"), useBytes = TRUE)
message("Cross-cohort meta-analysis and shared pathway analysis completed")
