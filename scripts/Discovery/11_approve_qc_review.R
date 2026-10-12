#!/usr/bin/env Rscript

root <- Sys.getenv("BD_ANALYSIS_ROOT", "inputs/discovery/bd_rhapsody_10/downstream/10samples_scrna")
result <- file.path(root, "results", "01_qc")
samples <- c("D01", "D02", "D03", "D04", "D05", "D06", "D07", "D08", "D09", "D10")
decision_id <- paste0("qc_exception_review_", format(Sys.time(), "%Y%m%dT%H%M%S%z"))
approved_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")

write_csv_atomic <- function(value, path, gzip = FALSE) {
  temporary <- paste0(path, ".tmp.", Sys.getpid())
  if (gzip) {
    connection <- gzfile(temporary, open = "wt")
    on.exit(close(connection), add = TRUE)
    write.csv(value, connection, row.names = FALSE)
    close(connection)
    on.exit(NULL, add = FALSE)
  } else {
    write.csv(value, temporary, row.names = FALSE)
  }
  if (!file.rename(temporary, path)) stop("Could not atomically publish ", path)
}

thresholds <- read.csv(file.path(result, "qc_thresholds_proposed.csv"), stringsAsFactors = FALSE)
summary <- read.csv(file.path(result, "10samples_qc_summary_pre_filter.csv"), stringsAsFactors = FALSE)
guardrails <- read.csv(file.path(result, "qc_guardrails.csv"), stringsAsFactors = FALSE)
cell_qc <- read.csv(gzfile(file.path(result, "cell_qc_pre_filter.csv.gz")), stringsAsFactors = FALSE)

if (!identical(thresholds$sample_id, samples)) stop("Proposed threshold sample order changed")
if (!identical(summary$sample_id, samples)) stop("QC summary sample order changed")
if (nrow(cell_qc) != sum(summary$cells_input)) stop("Cell-level QC row count changed")

failures <- guardrails[guardrails$status == "FAIL", , drop = FALSE]
expected_failure_keys <- sort(c(
  "D01|emptydrops_only_ratio", "D02|emptydrops_only_ratio",
  "D03|sample_robust_z", "D04|sample_robust_z", "D04|scrublet_doublet_rate",
  "D05|emptydrops_only_ratio", "D06|emptydrops_only_ratio",
  "D07|emptydrops_only_ratio", "D07|sample_robust_z",
  "D08|emptydrops_only_ratio", "D09|emptydrops_only_ratio", "D10|emptydrops_only_ratio"
))
actual_failure_keys <- sort(paste(failures$sample_id, failures$guardrail, sep = "|"))
if (!identical(actual_failure_keys, expected_failure_keys)) {
  stop("Guardrail failure set changed; refusing approved exception")
}

cell_qc$qc_keep_proposed <- as.logical(cell_qc$qc_keep_proposed)
cell_qc$predicted_doublet <- as.logical(cell_qc$predicted_doublet)
if (anyNA(cell_qc$qc_keep_proposed) || anyNA(cell_qc$predicted_doublet)) stop("Invalid logical QC fields")
if (any(cell_qc$qc_keep_proposed & cell_qc$predicted_doublet)) stop("Proposed keep set contains Scrublet doublets")

decision <- failures
decision$decision <- "APPROVED_EXCEPTION"
decision$action <- ifelse(
  decision$guardrail == "emptydrops_only_ratio",
  "retain_BD_filtered_only_do_not_add_emptyDrops_only",
  ifelse(
    decision$guardrail == "scrublet_doublet_rate",
    "remove_predicted_doublets_retain_patient",
    "retain_patient_add_quality_flag_and_sensitivity_analysis"
  )
)
decision$rationale <- ifelse(
  decision$guardrail == "emptydrops_only_ratio",
  "emptyDrops is diagnostic; canonical BD filtered RSEC MEX remains the cell universe",
  ifelse(
    decision$guardrail == "scrublet_doublet_rate",
    "all predicted doublets are excluded and the retained patient cell count remains above guardrails",
    "patient is retained as the statistical replicate and tested in prespecified sensitivity analyses"
  )
)
decision$decision_id <- decision_id
decision$approved_at <- approved_at

sample_flags <- data.frame(
  sample_id = samples,
  qc_outlier_flag = samples %in% c("D03", "D04", "D07"),
  high_doublet_flag = samples == "D04",
  primary_analysis = "INCLUDE",
  sensitivity_action = ifelse(samples %in% c("D03", "D04", "D07"), "PRESPECIFIED_EXCLUSION", "NONE"),
  stringsAsFactors = FALSE
)
sensitivity_plan <- data.frame(
  scenario_id = c("primary_all_10", "sensitivity_exclude_NS4", "sensitivity_exclude_all_qc_outliers"),
  excluded_samples = c("", "D04", "D03;D04;D07"),
  rationale = c(
    "Primary patient-replicate analysis retains all ten samples",
    "Tests influence of the high-doublet and strongest robust-z outlier",
    "Tests joint influence of all robust-z flagged samples"
  ),
  min_patients_per_severity = c(3L, 3L, 3L),
  required_stage = "04_compare",
  status = "PRESPECIFIED",
  stringsAsFactors = FALSE
)

approved <- thresholds
approved$approval_mode <- "EXPLICIT_EXCEPTION_REVIEW"
approved$approved_at <- approved_at
approved$decision_id <- decision_id

cell_qc$qc_keep <- cell_qc$qc_keep_proposed
retained <- table(cell_qc$sample_id[cell_qc$qc_keep])
summary$cells_qc_keep <- as.integer(retained[match(summary$sample_id, names(retained))])
summary$retention_rate <- summary$cells_qc_keep / summary$cells_input
summary$qc_exception_flag <- summary$sample_id %in% c("D03", "D04", "D07")
if (!all(summary$cells_qc_keep == summary$cells_proposed_keep)) stop("Reviewed keep counts differ from proposed counts")

approval <- data.frame(
  approval_mode = "EXPLICIT_EXCEPTION_REVIEW",
  status = "APPROVED",
  failed_guardrails = nrow(failures),
  timestamp = approved_at,
  decision_id = decision_id,
  reviewer = "user_confirmed_via_codex",
  conditions = "BD_filtered_only;remove_Scrublet_doublets;retain_all_patients;Stage04_sensitivity",
  stringsAsFactors = FALSE
)
keep_columns <- cell_qc[cell_qc$qc_keep, c("sample_id", "barcode", "cell_id"), drop = FALSE]

write_csv_atomic(decision, file.path(result, "qc_review_decision.csv"))
write_csv_atomic(sample_flags, file.path(result, "qc_sample_flags.csv"))
write_csv_atomic(sensitivity_plan, file.path(result, "qc_sensitivity_analysis_plan.csv"))
write_csv_atomic(approved, file.path(result, "qc_thresholds_approved.csv"))
write_csv_atomic(approved, file.path(result, "qc_thresholds.csv"))
write_csv_atomic(summary, file.path(result, "10samples_qc_summary.csv"))
write_csv_atomic(cell_qc, file.path(result, "cell_qc.csv.gz"), gzip = TRUE)
write_csv_atomic(keep_columns, file.path(result, "cells_keep.csv.gz"), gzip = TRUE)
write_csv_atomic(approval, file.path(result, "qc_approval.csv"))

report_path <- file.path(result, "qc_report.md")
report <- readLines(report_path, warn = FALSE)
review_section <- c(
  "", "## Reviewed exception approval", "",
  paste0("- Decision ID: `", decision_id, "`"),
  paste0("- Approved at: ", approved_at),
  "- All ten patient samples remain in the primary analysis.",
  "- Canonical BD filtered RSEC MEX remains the cell universe; emptyDrops-only barcodes are not added.",
  "- Predicted Scrublet doublets are excluded, including all D04 calls.",
  "- D03, D04, and D07 carry QC outlier flags.",
  "- Stage 04 must run the prespecified D04-only and all-QC-outlier exclusion sensitivity analyses."
)
writeLines(c(report, review_section), report_path)
writeLines(capture.output(sessionInfo()), file.path(result, "qc_review_r_sessionInfo.txt"))
message("Explicit QC exception approval published: ", decision_id, "; retained cells=", nrow(keep_columns))
