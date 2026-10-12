#!/usr/bin/env Rscript

root <- Sys.getenv("BD_ANALYSIS_ROOT", "inputs/discovery/bd_rhapsody_10/downstream/10samples_scrna")
result <- file.path(root, "results", "04_compare")
composition <- read.csv(file.path(result, "sample_cell_type_composition.csv"), stringsAsFactors = FALSE)
composition_tests <- read.csv(file.path(result, "cell_type_composition_comparison.csv"), stringsAsFactors = FALSE)
de <- read.csv(gzfile(file.path(result, "all_cell_types_pseudobulk_de.csv.gz")), stringsAsFactors = FALSE)
summary <- read.csv(file.path(result, "pseudobulk_de_summary.csv"), stringsAsFactors = FALSE)
sample_counts <- read.csv(file.path(result, "pseudobulk_sample_cell_counts.csv"), stringsAsFactors = FALSE)
sensitivity_plan <- read.csv(file.path(root, "results", "01_qc", "qc_sensitivity_analysis_plan.csv"), stringsAsFactors = FALSE)
sensitivity_de <- read.csv(gzfile(file.path(result, "all_cell_types_pseudobulk_sensitivity_de.csv.gz")), stringsAsFactors = FALSE)
sensitivity_summary <- read.csv(file.path(result, "pseudobulk_sensitivity_summary.csv"), stringsAsFactors = FALSE)
sensitivity_concordance <- read.csv(file.path(result, "pseudobulk_sensitivity_concordance.csv"), stringsAsFactors = FALSE)
checks <- list()
add <- function(check, passed, detail) checks[[length(checks) + 1L]] <<- data.frame(check = check, status = if (passed) "PASS" else "FAIL", detail = detail)
add("composition_has_10_patients", length(unique(composition$sample_id)) == 10, paste(unique(composition$sample_id), collapse = ","))
sums <- aggregate(proportion ~ sample_id, composition, sum)
add("composition_sums_to_one", all(abs(sums$proportion - 1) < 1e-8), paste(range(sums$proportion), collapse = "-"))
required_composition <- c("mean_difference_severe_minus_mild", "bootstrap_ci_low", "bootstrap_ci_high", "wilcoxon_p_value", "fdr")
add("composition_effect_ci_fdr", all(required_composition %in% colnames(composition_tests)) && all(is.finite(as.matrix(composition_tests[, required_composition]))), paste(required_composition, collapse = ","))
required_de <- c("cell_type_broad", "gene", "logFC", "ci95_low", "ci95_high", "PValue", "FDR", "n_mild_patients", "n_severe_patients")
add("de_effect_ci_fdr", all(required_de %in% colnames(de)) && all(is.finite(de$logFC)) && all(is.finite(de$ci95_low)) && all(is.finite(de$ci95_high)) && all(is.finite(de$FDR)), paste0("rows=", nrow(de)))
add("patient_replicate_guardrail", all(summary$n_mild_patients >= 3 & summary$n_severe_patients >= 3), paste0("cell_types=", nrow(summary)))
add("no_cell_level_de_columns", !any(c("cell_id", "barcode") %in% colnames(de)), "pseudobulk-only result schema")
add("sample_cell_counts_present", nrow(sample_counts) > 0 && length(unique(sample_counts$sample_id)) == 10, paste0("rows=", nrow(sample_counts)))
expected_scenarios <- c("sensitivity_exclude_NS4", "sensitivity_exclude_all_qc_outliers")
plan_rows <- sensitivity_plan[sensitivity_plan$scenario_id %in% expected_scenarios, , drop = FALSE]
add(
  "sensitivity_plan_fixed",
  identical(plan_rows$scenario_id, expected_scenarios) && identical(plan_rows$excluded_samples, c("D04", "D03;D04;D07")),
  paste(plan_rows$scenario_id, plan_rows$excluded_samples, collapse = ";")
)
summary_coverage <- table(sensitivity_summary$cell_type_broad, sensitivity_summary$scenario_id)
add("sensitivity_scenario_coverage", nrow(summary_coverage) == nrow(summary) && all(summary_coverage[, expected_scenarios, drop = FALSE] == 1), paste0("rows=", nrow(sensitivity_summary)))
completed_sensitivity <- sensitivity_summary[sensitivity_summary$status == "COMPLETED", , drop = FALSE]
add("sensitivity_patient_guardrail", nrow(completed_sensitivity) > 0 && all(completed_sensitivity$n_mild_patients >= 3 & completed_sensitivity$n_severe_patients >= 3), paste0("completed=", nrow(completed_sensitivity)))
required_sensitivity_de <- c("scenario_id", "excluded_samples", "cell_type_broad", "gene", "logFC", "ci95_low", "ci95_high", "PValue", "FDR", "n_mild_patients", "n_severe_patients")
add("sensitivity_de_effect_ci_fdr", nrow(sensitivity_de) > 0 && all(required_sensitivity_de %in% colnames(sensitivity_de)) && all(is.finite(sensitivity_de$logFC)) && all(is.finite(sensitivity_de$ci95_low)) && all(is.finite(sensitivity_de$ci95_high)) && all(is.finite(sensitivity_de$FDR)), paste0("rows=", nrow(sensitivity_de)))
add("sensitivity_concordance", nrow(sensitivity_concordance) == nrow(completed_sensitivity) && all(is.finite(sensitivity_concordance$logFC_spearman)) && all(is.finite(sensitivity_concordance$direction_agreement)), paste0("rows=", nrow(sensitivity_concordance)))
enrichment <- read.csv(gzfile(file.path(result, "all_cell_types_GO_BP_cameraPR.csv.gz")), stringsAsFactors = FALSE)
required_enrichment <- c("go_id", "term", "term_source", "Direction", "PValue", "FDR")
add("enrichment_complete", nrow(enrichment) > 0 && all(required_enrichment %in% colnames(enrichment)), paste0("rows=", nrow(enrichment)))
add(
  "enrichment_term_provenance",
  all(!is.na(enrichment$term) & nzchar(enrichment$term)) &&
    all(enrichment$term == enrichment$go_id) &&
    all(enrichment$term_source == "GO_ID_offline"),
  "term_source=GO_ID_offline"
)
required <- c("comparison_report.md", "pseudobulk_de_summary.csv", "cell_type_composition_comparison.csv", "plots/cell_type_composition_by_patient.pdf", "all_cell_types_pseudobulk_sensitivity_de.csv.gz", "pseudobulk_sensitivity_summary.csv", "pseudobulk_sensitivity_concordance.csv")
missing <- required[!file.exists(file.path(result, required)) | file.info(file.path(result, required))$size == 0]
add("required_outputs", !length(missing), paste(missing, collapse = ","))
validation <- do.call(rbind, checks)
write.csv(validation, file.path(result, "04_compare_validation.csv"), row.names = FALSE)
writeLines(c("# Stage 04 validation", "", paste0("- ", validation$status, ": ", validation$check, " (", validation$detail, ")")), file.path(result, "04_compare_validation.md"))
print(validation, row.names = FALSE)
if (any(validation$status != "PASS")) quit(status = 1L)
