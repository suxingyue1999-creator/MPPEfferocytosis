#!/usr/bin/env Rscript

lib <- Sys.getenv("BD_R_LIB", "inputs/discovery/bd_rhapsody_10/downstream/10samples_scrna/env/host-Rlib")
.libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages({
  library(Matrix)
  library(edgeR)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(ggplot2)
})

options(stringsAsFactors = FALSE, warn = 1)

# Public adaptation: severity is read from coded metadata, not from sample labels.
discovery_metadata_file <- Sys.getenv("MPP_DISCOVERY_METADATA", "inputs/discovery.csv")
discovery_metadata_public <- read.csv(discovery_metadata_file, stringsAsFactors = FALSE)
if (anyDuplicated(discovery_metadata_public$sample_id)) stop("Duplicate coded discovery sample")
discovery_severity <- setNames(discovery_metadata_public$severity_en, discovery_metadata_public$sample_id)

set.seed(20260727)

cross_root <- Sys.getenv(
  "BD_CROSS_ROOT",
  "inputs/discovery/bd_rhapsody_10/downstream/10_23_cross_cohort_v1"
)
discovery_root <- Sys.getenv(
  "BD_DISCOVERY_ROOT",
  "inputs/discovery/bd_rhapsody_10/downstream/10samples_scrna_mex_v4"
)
stage02 <- file.path(cross_root, "results", "02_harmonized_validation")
stage03 <- file.path(cross_root, "results", "03_meta_pathway")
result <- file.path(cross_root, "results", "04_replicated_mechanism")
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

replication <- read_csv_path(file.path(stage03, "lineage_replication_decisions.csv"))
eligible_lineages <- replication$harmonized_lineage[replication$mechanism_eligible %in% c(TRUE, "TRUE")]
if (!identical(eligible_lineages, "Macrophage_Monocyte")) {
  stop("Mechanism scope changed; expected only Macrophage_Monocyte but found: ", paste(eligible_lineages, collapse = ", "))
}
pathway_meta <- read_gz_csv_path(file.path(stage03, "two_cohort_pathway_integration.csv.gz"))
strict_pathways <- pathway_meta[
  pathway_meta$harmonized_lineage == "Macrophage_Monocyte" &
    pathway_meta$strictly_replicated_pathway %in% c(TRUE, "TRUE"),
  ,
  drop = FALSE
]
if (nrow(strict_pathways) < 1L) stop("No strictly replicated macrophage/monocyte pathway is available")
gene_meta <- read_gz_csv_path(file.path(stage03, "two_cohort_gene_meta_analysis.csv.gz"))
gene_meta <- gene_meta[gene_meta$harmonized_lineage == "Macrophage_Monocyte", , drop = FALSE]
common_genes <- sort(unique(gene_meta$gene))

symbol_to_entrez <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = common_genes,
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
pathway_members <- do.call(rbind, lapply(seq_len(nrow(strict_pathways)), function(index) {
  pathway <- strict_pathways[index, , drop = FALSE]
  entrez_ids <- unique(go_map$ENTREZID[go_map$GO == pathway$go_id])
  genes <- names(symbol_to_entrez)[!is.na(symbol_to_entrez) & symbol_to_entrez %in% entrez_ids]
  data.frame(
    go_id = pathway$go_id,
    term = pathway$term,
    pathway_direction = pathway$discovery_direction,
    gene = genes,
    stringsAsFactors = FALSE
  )
}))
pathway_members <- unique(pathway_members)
write_csv(pathway_members, "strict_pathway_gene_membership.csv")

membership_counts <- aggregate(gene ~ go_id + term, pathway_members, length)
names(membership_counts)[[3]] <- "reconstructed_NGenes"
membership_counts$cameraPR_NGenes <- strict_pathways$NGenes[match(membership_counts$go_id, strict_pathways$go_id)]
membership_counts$exact_match <- membership_counts$reconstructed_NGenes == membership_counts$cameraPR_NGenes
write_csv(membership_counts, "strict_pathway_membership_audit.csv")

count_root <- file.path(stage02, "pseudobulk_counts")
pb_inputs <- list(
  discovery_10 = list(
    counts = readRDS(file.path(count_root, "discovery_harmonized_raw_counts.rds")),
    meta = read_csv_path(file.path(count_root, "discovery_harmonized_metadata.csv"))
  ),
  validation_23 = list(
    counts = readRDS(file.path(count_root, "validation_harmonized_raw_counts.rds")),
    meta = read_csv_path(file.path(count_root, "validation_harmonized_metadata.csv"))
  )
)

score_frames <- list()
score_index <- 0L
for (cohort in names(pb_inputs)) {
  pb <- pb_inputs[[cohort]]
  pb$counts <- as(pb$counts, "dgCMatrix")
  columns <- which(pb$meta$harmonized_lineage == "Macrophage_Monocyte" & pb$meta$cell_count >= 50L)
  meta <- pb$meta[columns, , drop = FALSE]
  counts <- pb$counts[, columns, drop = FALSE]
  meta <- meta[match(colnames(counts), meta$pseudobulk_id), , drop = FALSE]
  meta$severity <- factor(as.character(meta$severity), levels = c("mild", "severe"))
  dge <- calcNormFactors(DGEList(counts = counts), method = "TMM")
  log_cpm <- cpm(dge, log = TRUE, prior.count = 1)

  for (index in seq_len(nrow(strict_pathways))) {
    pathway <- strict_pathways[index, , drop = FALSE]
    genes <- pathway_members$gene[pathway_members$go_id == pathway$go_id]
    genes <- intersect(genes, rownames(log_cpm))
    gene_scaled <- t(scale(t(log_cpm[genes, , drop = FALSE])))
    gene_scaled[!is.finite(gene_scaled)] <- 0
    raw_score <- colMeans(gene_scaled)
    standardized_score <- as.numeric(scale(raw_score))
    score_index <- score_index + 1L
    score_frames[[score_index]] <- data.frame(
      cohort = cohort,
      sample_id = meta$sample_id,
      severity = as.character(meta$severity),
      go_id = pathway$go_id,
      term = pathway$term,
      pathway_direction = pathway$discovery_direction,
      genes_scored = length(genes),
      raw_score = raw_score,
      standardized_score = standardized_score,
      batch_wave = if ("batch_wave" %in% names(meta)) as.character(meta$batch_wave) else NA_character_,
      allergy_group = if ("allergy_group" %in% names(meta)) as.character(meta$allergy_group) else NA_character_,
      plastic_group = if ("plastic_group" %in% names(meta)) as.character(meta$plastic_group) else NA_character_,
      stringsAsFactors = FALSE
    )
  }
}
pathway_scores <- do.call(rbind, score_frames)
write_csv(pathway_scores, "replicated_pathway_scores_by_sample.csv")

fit_score_model <- function(frame, cohort, model_id, adjusted = FALSE) {
  frame$severity <- factor(frame$severity, levels = c("mild", "severe"))
  if (adjusted) {
    frame$batch_wave <- factor(frame$batch_wave)
    frame$allergy_group <- factor(frame$allergy_group, levels = c("no", "yes"))
    frame$plastic_group <- factor(frame$plastic_group, levels = c("no", "yes"))
  }
  formula <- if (adjusted) {
    standardized_score ~ batch_wave + severity + allergy_group + plastic_group
  } else {
    standardized_score ~ severity
  }
  output <- do.call(rbind, lapply(split(frame, frame$go_id), function(current) {
    fit <- lm(formula, data = current)
    if (qr(model.matrix(fit))$rank != ncol(model.matrix(fit))) stop("Pathway score model is not full rank")
    coefficients <- summary(fit)$coefficients
    coefficients <- coefficients[rownames(coefficients) != "(Intercept)", , drop = FALSE]
    data.frame(
      cohort = cohort,
      model_id = model_id,
      go_id = current$go_id[[1]],
      term = current$term[[1]],
      pathway_direction = current$pathway_direction[[1]],
      coefficient = rownames(coefficients),
      estimate = coefficients[, "Estimate"],
      standard_error = coefficients[, "Std. Error"],
      statistic = coefficients[, "t value"],
      p_value = coefficients[, "Pr(>|t|)"],
      samples_used = nrow(current),
      stringsAsFactors = FALSE
    )
  }))
  output$fdr_by_coefficient <- ave(output$p_value, output$coefficient, FUN = function(x) p.adjust(x, "BH"))
  output
}

score_models <- rbind(
  fit_score_model(pathway_scores[pathway_scores$cohort == "discovery_10", ], "discovery_10", "severity_only"),
  fit_score_model(pathway_scores[pathway_scores$cohort == "validation_23", ], "validation_23", "severity_only"),
  fit_score_model(pathway_scores[pathway_scores$cohort == "validation_23", ], "validation_23", "covariate_adjusted", adjusted = TRUE)
)
write_csv(score_models, "replicated_pathway_score_models.csv")

meta_two <- function(effect1, se1, effect2, se2) {
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
  random_effect <- (random_weight1 * effect1 + random_weight2 * effect2) / (random_weight1 + random_weight2)
  random_se <- sqrt(1 / (random_weight1 + random_weight2))
  data.frame(
    fixed_effect = fixed_effect,
    fixed_se = fixed_se,
    fixed_p_value = 2 * pnorm(-abs(fixed_effect / fixed_se)),
    q = q,
    i2_percent = ifelse(q > 0, pmax(0, (q - 1) / q) * 100, 0),
    tau2 = tau2,
    random_effect = random_effect,
    random_se = random_se,
    random_p_value = 2 * pnorm(-abs(random_effect / random_se)),
    stringsAsFactors = FALSE
  )
}

discovery_score_effect <- score_models[
  score_models$cohort == "discovery_10" & score_models$model_id == "severity_only" & score_models$coefficient == "severitysevere",
  c("go_id", "term", "pathway_direction", "estimate", "standard_error", "p_value", "fdr_by_coefficient")
]
validation_score_effect <- score_models[
  score_models$cohort == "validation_23" & score_models$model_id == "severity_only" & score_models$coefficient == "severitysevere",
  c("go_id", "estimate", "standard_error", "p_value", "fdr_by_coefficient")
]
score_meta <- merge(discovery_score_effect, validation_score_effect, by = "go_id", suffixes = c("_discovery", "_validation"))
score_meta <- cbind(
  score_meta,
  meta_two(
    score_meta$estimate_discovery,
    score_meta$standard_error_discovery,
    score_meta$estimate_validation,
    score_meta$standard_error_validation
  )
)
score_meta$fixed_fdr <- p.adjust(score_meta$fixed_p_value, "BH")
score_meta$random_fdr <- p.adjust(score_meta$random_p_value, "BH")
score_meta$direction_concordant <- sign(score_meta$estimate_discovery) == sign(score_meta$estimate_validation)
write_csv(score_meta, "replicated_pathway_score_meta_analysis.csv")

pathway_gene_evidence <- merge(pathway_members, gene_meta, by = "gene", all.x = TRUE)
pathway_gene_evidence$direction_sign <- ifelse(pathway_gene_evidence$pathway_direction == "Up", 1, -1)
pathway_gene_evidence$direction_supported_both <- with(
  pathway_gene_evidence,
  sign(discovery_effect_log2fc) == direction_sign & sign(validation_effect_log2fc) == direction_sign
)
pathway_gene_evidence$aligned_meta_score <- pathway_gene_evidence$direction_sign * pathway_gene_evidence$fixed_z
pathway_gene_evidence <- pathway_gene_evidence[order(pathway_gene_evidence$go_id, -pathway_gene_evidence$aligned_meta_score), , drop = FALSE]
write.csv(pathway_gene_evidence, gzfile(file.path(result, "replicated_pathway_gene_evidence.csv.gz")), row.names = FALSE)
top_pathway_genes <- do.call(rbind, lapply(split(pathway_gene_evidence, pathway_gene_evidence$go_id), function(frame) head(frame, 10L)))
write_csv(top_pathway_genes, "top_replicated_pathway_genes.csv")

stage02_inventory <- read_csv_path(file.path(stage02, "stage_02_input_inventory.csv"))
validation_root <- sub(
  "/results/.*$",
  "",
  stage02_inventory$path[grepl("pseudobulk_raw_counts[.]rds$", stage02_inventory$path)]
)
validation_subtype_path <- file.path(
  validation_root,
  "results",
  "07_lineage_pseudobulk",
  "Macrophage_Monocyte",
  "subtype_composition_by_sample.csv"
)
validation_subtypes <- read_csv_path(validation_subtype_path)
validation_subtypes$severity <- ifelse(validation_subtypes$severity == "\u91cd\u75c7", "severe", "mild")
validation_subtypes$allergy_group <- ifelse(validation_subtypes$allergy_group == "\u662f", "yes", "no")
validation_subtypes$plastic_group <- ifelse(validation_subtypes$plastic_group == "\u662f", "yes", "no")

discovery_embedding_path <- file.path(
  discovery_root,
  "results",
  "03_annotate",
  "lineage_recluster_embeddings.csv.gz"
)
discovery_embeddings <- read_gz_csv_path(discovery_embedding_path)
discovery_embeddings <- discovery_embeddings[discovery_embeddings$lineage == "Macrophage_Monocyte", , drop = FALSE]
discovery_embeddings$sample_id <- sub("_[^_]+$", "", discovery_embeddings$cell_id)
discovery_counts <- aggregate(cell_id ~ sample_id + cell_type_subtype, discovery_embeddings, length)
names(discovery_counts)[[3]] <- "cell_count"
discovery_samples <- sort(unique(discovery_embeddings$sample_id))
discovery_labels <- sort(unique(discovery_embeddings$cell_type_subtype))
discovery_grid <- expand.grid(sample_id = discovery_samples, cell_type_subtype = discovery_labels, stringsAsFactors = FALSE)
discovery_subtypes <- merge(discovery_grid, discovery_counts, by = c("sample_id", "cell_type_subtype"), all.x = TRUE)
discovery_subtypes$cell_count[is.na(discovery_subtypes$cell_count)] <- 0L
discovery_totals <- aggregate(cell_count ~ sample_id, discovery_subtypes, sum)
names(discovery_totals)[[2]] <- "lineage_total"
discovery_subtypes <- merge(discovery_subtypes, discovery_totals, by = "sample_id")
discovery_subtypes$severity <- unname(discovery_severity[discovery_subtypes$sample_id])
discovery_subtypes$fraction <- discovery_subtypes$cell_count / discovery_subtypes$lineage_total

subtype_mapping <- data.frame(
  discovery_label = c("Alveolar_Macrophage", "Cycling_Myeloid"),
  validation_label = c("FABP4_alveolar", "proliferating"),
  harmonized_state = c("Alveolar_macrophage", "Cycling_myeloid"),
  meta_eligible = TRUE,
  mapping_basis = c(
    "FABP4/alveolar macrophage identity",
    "cell-cycle/proliferation identity"
  ),
  stringsAsFactors = FALSE
)
write_csv(subtype_mapping, "macrophage_subtype_mapping.csv")

make_state_frame <- function(frame, cohort, label_column, mapping_column) {
  output <- do.call(rbind, lapply(seq_len(nrow(subtype_mapping)), function(index) {
    label <- subtype_mapping[[mapping_column]][[index]]
    current <- frame[frame[[label_column]] == label, , drop = FALSE]
    data.frame(
      cohort = cohort,
      sample_id = current$sample_id,
      severity = current$severity,
      harmonized_state = subtype_mapping$harmonized_state[[index]],
      state_count = current$cell_count,
      lineage_total = current$lineage_total,
      fraction = current$fraction,
      batch_wave = if ("batch_wave" %in% names(current)) as.character(current$batch_wave) else NA_character_,
      allergy_group = if ("allergy_group" %in% names(current)) as.character(current$allergy_group) else NA_character_,
      plastic_group = if ("plastic_group" %in% names(current)) as.character(current$plastic_group) else NA_character_,
      stringsAsFactors = FALSE
    )
  }))
  output
}

state_composition <- rbind(
  make_state_frame(discovery_subtypes, "discovery_10", "cell_type_subtype", "discovery_label"),
  make_state_frame(validation_subtypes, "validation_23", "subtype", "validation_label")
)
state_composition <- state_composition[
  !(state_composition$cohort == "discovery_10" & state_composition$lineage_total < 50L),
  ,
  drop = FALSE
]
write_csv(state_composition, "harmonized_macrophage_state_composition.csv")

fit_state_model <- function(frame, cohort, model_id, adjusted = FALSE) {
  frame$severity <- factor(frame$severity, levels = c("mild", "severe"))
  if (adjusted) {
    frame$batch_wave <- factor(frame$batch_wave)
    frame$allergy_group <- factor(frame$allergy_group, levels = c("no", "yes"))
    frame$plastic_group <- factor(frame$plastic_group, levels = c("no", "yes"))
  }
  output <- do.call(rbind, lapply(split(frame, frame$harmonized_state), function(current) {
    formula <- if (adjusted) {
      cbind(state_count, lineage_total - state_count) ~ batch_wave + severity + allergy_group + plastic_group
    } else {
      cbind(state_count, lineage_total - state_count) ~ severity
    }
    fit <- glm(formula, data = current, family = quasibinomial())
    if (qr(model.matrix(fit))$rank != ncol(model.matrix(fit))) stop("State-composition model is not full rank")
    coefficients <- summary(fit)$coefficients
    coefficients <- coefficients[rownames(coefficients) != "(Intercept)", , drop = FALSE]
    data.frame(
      cohort = cohort,
      model_id = model_id,
      harmonized_state = current$harmonized_state[[1]],
      coefficient = rownames(coefficients),
      estimate_log_odds = coefficients[, "Estimate"],
      standard_error = coefficients[, "Std. Error"],
      statistic = coefficients[, "t value"],
      p_value = coefficients[, "Pr(>|t|)"],
      samples_used = nrow(current),
      stringsAsFactors = FALSE
    )
  }))
  output$fdr_by_coefficient <- ave(output$p_value, output$coefficient, FUN = function(x) p.adjust(x, "BH"))
  output
}

state_models <- rbind(
  fit_state_model(state_composition[state_composition$cohort == "discovery_10", ], "discovery_10", "severity_only"),
  fit_state_model(state_composition[state_composition$cohort == "validation_23", ], "validation_23", "severity_only"),
  fit_state_model(state_composition[state_composition$cohort == "validation_23", ], "validation_23", "covariate_adjusted", adjusted = TRUE)
)
write_csv(state_models, "harmonized_macrophage_state_models.csv")

state_d <- state_models[state_models$cohort == "discovery_10" & state_models$model_id == "severity_only" & state_models$coefficient == "severitysevere", ]
state_v <- state_models[state_models$cohort == "validation_23" & state_models$model_id == "severity_only" & state_models$coefficient == "severitysevere", ]
state_meta <- merge(
  state_d[, c("harmonized_state", "estimate_log_odds", "standard_error", "p_value", "fdr_by_coefficient")],
  state_v[, c("harmonized_state", "estimate_log_odds", "standard_error", "p_value", "fdr_by_coefficient")],
  by = "harmonized_state",
  suffixes = c("_discovery", "_validation")
)
state_meta <- cbind(
  state_meta,
  meta_two(
    state_meta$estimate_log_odds_discovery,
    state_meta$standard_error_discovery,
    state_meta$estimate_log_odds_validation,
    state_meta$standard_error_validation
  )
)
state_meta$fixed_fdr <- p.adjust(state_meta$fixed_p_value, "BH")
state_meta$direction_concordant <- sign(state_meta$estimate_log_odds_discovery) == sign(state_meta$estimate_log_odds_validation)
write_csv(state_meta, "harmonized_macrophage_state_meta_analysis.csv")

scope <- replication[, c("harmonized_lineage", "mechanism_eligible", "replication_status"), drop = FALSE]
scope$mechanism_executed <- scope$mechanism_eligible
scope$analysis_scope <- ifelse(
  scope$mechanism_eligible %in% c(TRUE, "TRUE"),
  "replicated_pathway_scores;clinical_covariate_sensitivity;harmonized_state_composition",
  "none"
)
scope$reason <- ifelse(
  scope$mechanism_eligible %in% c(TRUE, "TRUE"),
  "lineage passed the prespecified strong replication rule",
  "lineage did not pass the prespecified strong replication rule"
)
write_csv(scope, "mechanism_scope_by_lineage.csv")

method_decisions <- data.frame(
  method = c("pathway_state_scores", "clinical_covariate_sensitivity", "subtype_composition", "regulatory_network", "trajectory", "CellChat", "ADT_support"),
  decision = c("RUN", "RUN", "RUN", "SKIP", "SKIP", "SKIP", "SKIP"),
  reason = c(
    "seven pathways replicated with FDR < 0.05 in both cohorts",
    "validation cohort has a full-rank batch/allergy/plastic-adjusted design",
    "two macrophage states have explicit cross-cohort label mappings",
    "no prespecified cross-cohort regulon resource and no macrophage gene-level replication",
    "cross-sectional independent cohorts do not establish a shared temporal trajectory",
    "historical network is pooled and no sample-level cross-cohort sender-receiver test is available",
    "ADT is available only in validation and cannot be independently confirmed in discovery"
  ),
  stringsAsFactors = FALSE
)
write_csv(method_decisions, "mechanism_method_decisions.csv")

score_plot <- ggplot(pathway_scores, aes(severity, standardized_score, color = severity)) +
  geom_boxplot(outlier.shape = NA, color = "#555555", width = 0.6) +
  geom_jitter(width = 0.12, height = 0, alpha = 0.75, size = 1.2) +
  facet_grid(cohort ~ term, scales = "free_y") +
  scale_color_manual(values = c(mild = "#2878B5", severe = "#D64B40")) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  labs(x = NULL, y = "Standardized pathway score")
ggsave(file.path(figure_root, "replicated_pathway_scores.pdf"), score_plot, width = 18, height = 7, useDingbats = FALSE)
ggsave(file.path(figure_root, "replicated_pathway_scores.png"), score_plot, width = 18, height = 7, dpi = 180)

score_forest <- rbind(
  data.frame(go_id = score_meta$go_id, term = score_meta$term, estimate = "Discovery", effect = score_meta$estimate_discovery, se = score_meta$standard_error_discovery),
  data.frame(go_id = score_meta$go_id, term = score_meta$term, estimate = "Validation", effect = score_meta$estimate_validation, se = score_meta$standard_error_validation),
  data.frame(go_id = score_meta$go_id, term = score_meta$term, estimate = "Fixed meta", effect = score_meta$fixed_effect, se = score_meta$fixed_se)
)
score_forest$low <- score_forest$effect - 1.96 * score_forest$se
score_forest$high <- score_forest$effect + 1.96 * score_forest$se
score_forest$term <- factor(score_forest$term, levels = rev(unique(score_meta$term)))
score_forest_plot <- ggplot(score_forest, aes(term, effect, color = estimate)) +
  geom_hline(yintercept = 0, color = "#777777", linewidth = 0.3) +
  geom_errorbar(aes(ymin = low, ymax = high), position = position_dodge(width = 0.55), width = 0.18) +
  geom_point(position = position_dodge(width = 0.55), size = 1.6) +
  coord_flip() +
  scale_color_manual(values = c("Discovery" = "#2878B5", "Validation" = "#D64B40", "Fixed meta" = "#222222")) +
  theme_bw(base_size = 9) +
  labs(x = NULL, y = "Standardized severe-vs-mild effect", color = NULL)
ggsave(file.path(figure_root, "replicated_pathway_score_forest.pdf"), score_forest_plot, width = 11, height = 6, useDingbats = FALSE)
ggsave(file.path(figure_root, "replicated_pathway_score_forest.png"), score_forest_plot, width = 11, height = 6, dpi = 180)

state_plot <- ggplot(state_composition, aes(severity, fraction, color = severity)) +
  geom_boxplot(outlier.shape = NA, color = "#555555", width = 0.6) +
  geom_jitter(width = 0.12, height = 0, alpha = 0.75, size = 1.3) +
  facet_grid(cohort ~ harmonized_state, scales = "free_y") +
  scale_color_manual(values = c(mild = "#2878B5", severe = "#D64B40")) +
  theme_bw(base_size = 10) +
  theme(legend.position = "none") +
  labs(x = NULL, y = "Fraction within macrophage/monocyte lineage")
ggsave(file.path(figure_root, "harmonized_macrophage_state_composition.pdf"), state_plot, width = 10, height = 6, useDingbats = FALSE)
ggsave(file.path(figure_root, "harmonized_macrophage_state_composition.png"), state_plot, width = 10, height = 6, dpi = 180)

gene_plot_frame <- rbind(
  data.frame(go_id = top_pathway_genes$go_id, term = top_pathway_genes$term, gene = top_pathway_genes$gene, cohort = "Discovery", effect = top_pathway_genes$discovery_effect_log2fc),
  data.frame(go_id = top_pathway_genes$go_id, term = top_pathway_genes$term, gene = top_pathway_genes$gene, cohort = "Validation", effect = top_pathway_genes$validation_effect_log2fc),
  data.frame(go_id = top_pathway_genes$go_id, term = top_pathway_genes$term, gene = top_pathway_genes$gene, cohort = "Fixed meta", effect = top_pathway_genes$fixed_effect)
)
gene_effect_plot <- ggplot(gene_plot_frame, aes(cohort, gene, fill = effect)) +
  geom_tile(color = "white", linewidth = 0.15) +
  facet_wrap(~ term, scales = "free_y") +
  scale_fill_gradient2(low = "#2878B5", mid = "white", high = "#C43C39", midpoint = 0) +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = NULL, fill = "log2 effect")
ggsave(file.path(figure_root, "replicated_pathway_gene_effects.pdf"), gene_effect_plot, width = 14, height = 10, useDingbats = FALSE)
ggsave(file.path(figure_root, "replicated_pathway_gene_effects.png"), gene_effect_plot, width = 14, height = 10, dpi = 180)

report <- c(
  "# Replicated-lineage mechanism analysis",
  "",
  paste0("- Generated: ", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
  "- The mechanism scope is restricted to Macrophage_Monocyte because it is the only strongly replicated lineage.",
  paste0("- Strictly replicated GO-BP pathways analyzed: ", nrow(strict_pathways), "."),
  "- Pathway scores are patient-level standardized summaries derived from TMM log-CPM pseudobulks.",
  "- Validation score models include the prespecified batch wave, allergy, and plastic-surgery covariates.",
  "- Pathway score tests are supportive characterizations of already selected replicated pathways, not independent discovery tests.",
  "- Cross-cohort subtype composition is limited to alveolar macrophage and cycling myeloid states with explicit label mappings.",
  "- Regulatory-network, trajectory, CellChat, and ADT claims are not run because their cross-cohort inputs or statistical replicate structure are insufficient.",
  "- All state-composition models use patient samples as replicates; no cell is treated as an independent replicate."
)
writeLines(report, file.path(result, "replicated_mechanism_report.md"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(result, "sessionInfo.txt"), useBytes = TRUE)
message("Replicated-lineage mechanism analysis completed")
