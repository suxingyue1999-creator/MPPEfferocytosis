#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

cross_root <- Sys.getenv(
  "BD_CROSS_ROOT",
  "inputs/discovery/bd_rhapsody_10/downstream/10_23_cross_cohort_v1"
)
result <- file.path(cross_root, "results", "04_replicated_mechanism")

read_csv <- function(path) {
  if (!file.exists(path)) stop("Missing Stage 04 output: ", path)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8-BOM")
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

scope <- read_csv(file.path(result, "mechanism_scope_by_lineage.csv"))
decisions <- read_csv(file.path(result, "mechanism_method_decisions.csv"))
membership <- read_csv(file.path(result, "strict_pathway_gene_membership.csv"))
membership_audit <- read_csv(file.path(result, "strict_pathway_membership_audit.csv"))
scores <- read_csv(file.path(result, "replicated_pathway_scores_by_sample.csv"))
score_models <- read_csv(file.path(result, "replicated_pathway_score_models.csv"))
score_meta <- read_csv(file.path(result, "replicated_pathway_score_meta_analysis.csv"))
subtype_mapping <- read_csv(file.path(result, "macrophage_subtype_mapping.csv"))
state_composition <- read_csv(file.path(result, "harmonized_macrophage_state_composition.csv"))
state_models <- read_csv(file.path(result, "harmonized_macrophage_state_models.csv"))
state_meta <- read_csv(file.path(result, "harmonized_macrophage_state_meta_analysis.csv"))

executed <- scope$harmonized_lineage[scope$mechanism_executed %in% c(TRUE, "TRUE")]
add_check("mechanism_scope_only_replicated_lineage", identical(executed, "Macrophage_Monocyte"), paste(executed, collapse = ";"))
add_check(
  "strict_pathway_membership_exact",
  nrow(membership_audit) == 7L && all(membership_audit$exact_match %in% c(TRUE, "TRUE")),
  paste0("pathways=", nrow(membership_audit), ";genes=", nrow(membership))
)
score_keys <- unique(scores[, c("cohort", "go_id")])
add_check(
  "pathway_scores_complete",
  nrow(score_keys) == 14L && nrow(scores) == 224L && all(is.finite(scores$standardized_score)),
  paste0("rows=", nrow(scores), ";models=", nrow(score_keys))
)
primary_score_models <- score_models[score_models$model_id == "severity_only" & score_models$coefficient == "severitysevere", ]
adjusted_score_models <- score_models[score_models$model_id == "covariate_adjusted" & score_models$coefficient == "severitysevere", ]
add_check(
  "pathway_score_models_complete",
  nrow(primary_score_models) == 14L && nrow(adjusted_score_models) == 7L && all(score_models$standard_error > 0),
  paste0("primary=", nrow(primary_score_models), ";adjusted=", nrow(adjusted_score_models))
)
add_check(
  "pathway_score_meta_complete",
  nrow(score_meta) == 7L && all(is.finite(score_meta$fixed_effect) & score_meta$fixed_se > 0),
  paste0("pathways=", nrow(score_meta))
)
add_check(
  "subtype_mapping_prespecified",
  nrow(subtype_mapping) == 2L && all(subtype_mapping$meta_eligible %in% c(TRUE, "TRUE")) &&
    setequal(subtype_mapping$harmonized_state, c("Alveolar_macrophage", "Cycling_myeloid")),
  paste(subtype_mapping$harmonized_state, collapse = ";")
)
state_keys <- unique(state_composition[, c("cohort", "sample_id", "harmonized_state")])
add_check(
  "state_composition_patient_level",
  nrow(state_keys) == nrow(state_composition) && all(state_composition$fraction >= 0 & state_composition$fraction <= 1),
  paste0("rows=", nrow(state_composition))
)
primary_state_models <- state_models[state_models$model_id == "severity_only" & state_models$coefficient == "severitysevere", ]
adjusted_state_models <- state_models[state_models$model_id == "covariate_adjusted" & state_models$coefficient == "severitysevere", ]
add_check(
  "state_models_complete",
  nrow(primary_state_models) == 4L && nrow(adjusted_state_models) == 2L && all(state_models$standard_error > 0),
  paste0("primary=", nrow(primary_state_models), ";adjusted=", nrow(adjusted_state_models))
)
add_check(
  "state_meta_complete",
  nrow(state_meta) == 2L && all(is.finite(state_meta$fixed_effect) & state_meta$fixed_se > 0),
  paste0("states=", nrow(state_meta))
)
add_check(
  "unsupported_methods_skipped",
  all(decisions$decision[decisions$method %in% c("regulatory_network", "trajectory", "CellChat", "ADT_support")] == "SKIP"),
  paste(decisions$method[decisions$decision == "SKIP"], collapse = ";")
)

figures <- file.path(
  result,
  "figures",
  c(
    "replicated_pathway_scores.pdf", "replicated_pathway_scores.png",
    "replicated_pathway_score_forest.pdf", "replicated_pathway_score_forest.png",
    "harmonized_macrophage_state_composition.pdf", "harmonized_macrophage_state_composition.png",
    "replicated_pathway_gene_effects.pdf", "replicated_pathway_gene_effects.png"
  )
)
figure_info <- file.info(figures)
add_check("mechanism_figures_nonempty", all(file.exists(figures)) && all(figure_info$size > 5000), paste0("figures=", length(figures)))
required <- c("replicated_pathway_gene_evidence.csv.gz", "replicated_mechanism_report.md", "sessionInfo.txt")
missing <- required[!file.exists(file.path(result, required))]
add_check("required_stage_outputs", !length(missing), paste(missing, collapse = ";"))

validation <- do.call(rbind, checks)
write.csv(validation, file.path(result, "04_replicated_mechanism_validation.csv"), row.names = FALSE)
writeLines(
  c("# Stage 04 replicated-mechanism validation", "", capture.output(print(validation, row.names = FALSE))),
  file.path(result, "04_replicated_mechanism_validation.md"),
  useBytes = TRUE
)
if (any(validation$status != "PASS")) {
  stop("Stage 04 validation failed: ", paste(validation$check[validation$status != "PASS"], collapse = ", "))
}
message("Stage 04 replicated-mechanism validation PASS")
