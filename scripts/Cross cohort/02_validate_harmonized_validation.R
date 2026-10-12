#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

cross_root <- Sys.getenv(
  "BD_CROSS_ROOT",
  "inputs/discovery/bd_rhapsody_10/downstream/10_23_cross_cohort_v1"
)
result <- file.path(cross_root, "results", "02_harmonized_validation")

read_csv <- function(name) {
  path <- file.path(result, name)
  if (!file.exists(path)) stop("Missing Stage 02 output: ", path)
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

primary_summary <- read_csv("harmonized_primary_de_summary.csv")
adjusted_summary <- read_csv("validation_adjusted_de_summary.csv")
composition <- read_csv("harmonized_composition_by_sample.csv")
composition_tests <- read_csv("harmonized_composition_comparison.csv")
t_exclusion <- read_csv("discovery_T_lineage_exclusion_audit.csv")
t_support <- read_csv("validation_T_lineage_marker_support.csv")
mac_audit <- read_csv("validation_macrophage_monocyte_sum_audit.csv")
input_inventory <- read_csv("stage_02_input_inventory.csv")

expected_lineages <- c("B_lineage", "Macrophage_Monocyte", "Neutrophil", "T_lineage")
expected_primary <- expand.grid(
  cohort = c("discovery_10", "validation_23"),
  harmonized_lineage = expected_lineages,
  stringsAsFactors = FALSE
)
observed_primary <- unique(primary_summary[, c("cohort", "harmonized_lineage")])
primary_keys <- paste(observed_primary$cohort, observed_primary$harmonized_lineage)
expected_keys <- paste(expected_primary$cohort, expected_primary$harmonized_lineage)
add_check("four_lineages_two_cohorts", nrow(primary_summary) == 8L && setequal(primary_keys, expected_keys), paste0("models=", nrow(primary_summary)))
add_check("primary_models_completed", all(primary_summary$status == "COMPLETED"), paste(table(primary_summary$status), collapse = ";"))
add_check(
  "patient_replicate_guardrail",
  all(primary_summary$n_mild_patients >= 3L & primary_summary$n_severe_patients >= 3L & primary_summary$min_cells_per_pseudobulk >= 50L),
  paste0("min_group=", min(primary_summary$n_mild_patients, primary_summary$n_severe_patients), ";min_cells=", min(primary_summary$min_cells_per_pseudobulk))
)
add_check("primary_design_full_rank", all(primary_summary$design_rank == primary_summary$design_columns), paste0("models=", nrow(primary_summary)))
add_check(
  "effect_standard_errors_complete",
  all(primary_summary$finite_positive_effect_se == primary_summary$genes_tested),
  paste0("genes=", sum(primary_summary$genes_tested), ";se=", sum(primary_summary$finite_positive_effect_se))
)
add_check(
  "validation_adjusted_models_complete",
  nrow(adjusted_summary) == 4L && setequal(adjusted_summary$harmonized_lineage, expected_lineages) && all(adjusted_summary$status == "COMPLETED"),
  paste0("models=", nrow(adjusted_summary))
)
add_check("validation_adjusted_full_rank", all(adjusted_summary$design_rank == adjusted_summary$design_columns), paste0("models=", nrow(adjusted_summary)))
add_check(
  "macrophage_monocyte_summed_before_modeling",
  nrow(mac_audit) == 23L && all(mac_audit$exact_sum %in% c(TRUE, "TRUE")),
  paste0("samples=", nrow(mac_audit), ";exact=", sum(mac_audit$exact_sum %in% c(TRUE, "TRUE")))
)
add_check(
  "discovery_T_exclusions_recorded",
  nrow(t_exclusion) == 10L && all(t_exclusion$excluded_NK_Cycling_Unresolved >= 0L) && sum(t_exclusion$excluded_NK_Cycling_Unresolved) > 0L,
  paste0("samples=", nrow(t_exclusion), ";excluded=", sum(t_exclusion$excluded_NK_Cycling_Unresolved))
)
add_check(
  "validation_T_marker_support",
  nrow(t_support) == 1L && t_support$support_status == "PASS" && t_support$t_pseudobulks == 23L,
  paste0("status=", t_support$support_status, ";delta=", signif(t_support$median_score_delta, 3))
)
composition_keys <- paste(composition_tests$cohort, composition_tests$harmonized_lineage)
add_check("composition_models_complete", nrow(composition_tests) == 8L && setequal(composition_keys, expected_keys), paste0("tests=", nrow(composition_tests)))
add_check(
  "composition_patient_rows_complete",
  nrow(composition) == 132L && all(composition$proportion >= 0 & composition$proportion <= 1),
  paste0("rows=", nrow(composition))
)
add_check("required_inputs_present", nrow(input_inventory) == 4L && all(input_inventory$exists %in% c(TRUE, "TRUE")), paste0("files=", nrow(input_inventory)))

required_outputs <- c(
  "harmonized_primary_de.csv.gz",
  "validation_adjusted_de.csv.gz",
  "harmonized_validation_report.md",
  "sessionInfo.txt",
  file.path("pseudobulk_counts", "discovery_harmonized_raw_counts.rds"),
  file.path("pseudobulk_counts", "validation_harmonized_raw_counts.rds")
)
missing_outputs <- required_outputs[!file.exists(file.path(result, required_outputs))]
add_check("required_stage_outputs", !length(missing_outputs), paste(missing_outputs, collapse = ";"))

validation <- do.call(rbind, checks)
write.csv(validation, file.path(result, "02_harmonized_validation_checks.csv"), row.names = FALSE)
report <- c("# Stage 02 harmonized validation", "", capture.output(print(validation, row.names = FALSE)))
writeLines(report, file.path(result, "02_harmonized_validation_checks.md"), useBytes = TRUE)
if (any(validation$status != "PASS")) {
  stop("Stage 02 validation failed: ", paste(validation$check[validation$status != "PASS"], collapse = ", "))
}
message("Stage 02 harmonized validation PASS")
