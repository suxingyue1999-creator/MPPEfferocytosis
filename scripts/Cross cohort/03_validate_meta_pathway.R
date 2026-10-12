#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

cross_root <- Sys.getenv(
  "BD_CROSS_ROOT",
  "inputs/discovery/bd_rhapsody_10/downstream/10_23_cross_cohort_v1"
)
result <- file.path(cross_root, "results", "03_meta_pathway")

read_csv <- function(path) {
  if (!file.exists(path)) stop("Missing Stage 03 output: ", path)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8-BOM")
}

read_gz_csv <- function(path) {
  if (!file.exists(path)) stop("Missing Stage 03 output: ", path)
  read.csv(gzfile(path), stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8-BOM")
}

checks <- list()
add_check <- function(check, status, detail) {
  checks[[length(checks) + 1L]] <<- data.frame(
    check = check,
    status = if (isTRUE(status)) "PASS" else "FAIL",
    detail = as.character(detail),
    stringsAsFactors = FALSE
  )
}

gene_meta <- read_gz_csv(file.path(result, "two_cohort_gene_meta_analysis.csv.gz"))
concordance <- read_csv(file.path(result, "gene_effect_concordance_summary.csv"))
adjustment <- read_csv(file.path(result, "validation_adjustment_concordance.csv"))
pathway <- read_gz_csv(file.path(result, "shared_GO_BP_cameraPR.csv.gz"))
pathway_meta <- read_gz_csv(file.path(result, "two_cohort_pathway_integration.csv.gz"))
replication <- read_csv(file.path(result, "lineage_replication_decisions.csv"))
composition <- read_csv(file.path(result, "two_cohort_composition_meta_analysis.csv"))
expected_lineages <- c("B_lineage", "Macrophage_Monocyte", "Neutrophil", "T_lineage")

gene_counts <- table(gene_meta$harmonized_lineage)
add_check(
  "gene_meta_four_lineages",
  setequal(names(gene_counts), expected_lineages) && all(gene_counts >= 4000L),
  paste(names(gene_counts), gene_counts, collapse = ";")
)
add_check(
  "gene_meta_effects_complete",
  all(is.finite(gene_meta$fixed_effect) & gene_meta$fixed_se > 0 & is.finite(gene_meta$random_effect) & gene_meta$random_se > 0),
  paste0("rows=", nrow(gene_meta))
)
add_check(
  "gene_meta_probabilities_valid",
  all(gene_meta$fixed_p_value >= 0 & gene_meta$fixed_p_value <= 1 & gene_meta$fixed_fdr >= 0 & gene_meta$fixed_fdr <= 1 & gene_meta$q_p_value >= 0 & gene_meta$q_p_value <= 1),
  paste0("rows=", nrow(gene_meta))
)
add_check(
  "gene_replication_rule_enforced",
  all(!gene_meta$gene_replicated | (
    gene_meta$direction_concordant & gene_meta$discovery_FDR < 0.05 & gene_meta$validation_PValue < 0.05 &
      gene_meta$fixed_fdr < 0.05 & gene_meta$i2_percent < 75
  )),
  paste0("replicated=", sum(gene_meta$gene_replicated))
)
add_check("effect_concordance_complete", nrow(concordance) == 4L && setequal(concordance$harmonized_lineage, expected_lineages), paste0("lineages=", nrow(concordance)))
add_check("adjustment_concordance_complete", nrow(adjustment) == 4L && setequal(adjustment$harmonized_lineage, expected_lineages), paste0("lineages=", nrow(adjustment)))

pathway_keys <- unique(pathway[, c("cohort", "harmonized_lineage")])
add_check(
  "shared_pathway_models_complete",
  nrow(pathway_keys) == 8L && setequal(pathway_keys$cohort, c("discovery_10", "validation_23")) && setequal(pathway_keys$harmonized_lineage, expected_lineages),
  paste0("models=", nrow(pathway_keys))
)
universe_check <- aggregate(common_gene_universe ~ harmonized_lineage, pathway, function(x) length(unique(x)))
add_check("pathway_common_universe_locked", nrow(universe_check) == 4L && all(universe_check$common_gene_universe == 1L), paste0("lineages=", nrow(universe_check)))
add_check(
  "pathway_offline_method_locked",
  all(pathway$term_source == "GO_BP_org.Hs.eg.db_offline_common_universe"),
  paste(unique(pathway$term_source), collapse = ";")
)
add_check(
  "pathway_integration_complete",
  setequal(unique(pathway_meta$harmonized_lineage), expected_lineages) && all(pathway_meta$combined_p_value >= 0 & pathway_meta$combined_p_value <= 1),
  paste0("rows=", nrow(pathway_meta))
)
add_check(
  "strict_pathway_rule_enforced",
  all(!pathway_meta$strictly_replicated_pathway | (
    pathway_meta$direction_concordant & pathway_meta$discovery_fdr < 0.05 & pathway_meta$validation_fdr < 0.05
  )),
  paste0("replicated=", sum(pathway_meta$strictly_replicated_pathway))
)
add_check(
  "lineage_replication_decisions_complete",
  nrow(replication) == 4L && setequal(replication$harmonized_lineage, expected_lineages) && all(replication$replication_status %in% c("not_replicated", "replicated_limited", "replicated_strong")),
  paste(replication$harmonized_lineage, replication$replication_status, collapse = ";")
)
add_check(
  "mechanism_eligibility_rule_enforced",
  all(replication$mechanism_eligible %in% c(TRUE, "TRUE") == (replication$replicated_genes >= 3L | replication$strictly_replicated_pathways >= 1L)),
  paste(replication$harmonized_lineage[replication$mechanism_eligible %in% c(TRUE, "TRUE")], collapse = ";")
)
add_check("composition_meta_complete", nrow(composition) == 4L && setequal(composition$harmonized_lineage, expected_lineages), paste0("lineages=", nrow(composition)))

figures <- file.path(
  result,
  "figures",
  c(
    "gene_effect_concordance.pdf", "gene_effect_concordance.png",
    "top_gene_meta_forest.pdf", "top_gene_meta_forest.png",
    "shared_pathway_integration.pdf", "shared_pathway_integration.png",
    "composition_meta_forest.pdf", "composition_meta_forest.png"
  )
)
figure_info <- file.info(figures)
add_check("publication_figures_nonempty", all(file.exists(figures)) && all(figure_info$size > 5000), paste0("figures=", length(figures)))
required <- c("meta_pathway_report.md", "sessionInfo.txt")
missing <- required[!file.exists(file.path(result, required))]
add_check("required_stage_outputs", !length(missing), paste(missing, collapse = ";"))

validation <- do.call(rbind, checks)
write.csv(validation, file.path(result, "03_meta_pathway_validation.csv"), row.names = FALSE)
writeLines(
  c("# Stage 03 meta/pathway validation", "", capture.output(print(validation, row.names = FALSE))),
  file.path(result, "03_meta_pathway_validation.md"),
  useBytes = TRUE
)
if (any(validation$status != "PASS")) {
  stop("Stage 03 validation failed: ", paste(validation$check[validation$status != "PASS"], collapse = ", "))
}
message("Stage 03 meta/pathway validation PASS")
