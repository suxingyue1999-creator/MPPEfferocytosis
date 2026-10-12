#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

discovery_root <- Sys.getenv(
  "BD_DISCOVERY_ROOT",
  "inputs/discovery/bd_rhapsody_10/downstream/10samples_scrna_mex_v4"
)
cross_root <- Sys.getenv(
  "BD_CROSS_ROOT",
  "inputs/discovery/bd_rhapsody_10/downstream/10_23_cross_cohort_v1"
)
audit_root <- file.path(cross_root, "results", "00_audit")
result <- file.path(cross_root, "results", "01_discovery_freeze")
dir.create(result, recursive = TRUE, showWarnings = FALSE)

read_csv <- function(path) {
  if (!file.exists(path)) stop("Missing required input: ", path)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8-BOM")
}

write_csv <- function(frame, name) {
  write.csv(frame, file.path(result, name), row.names = FALSE, fileEncoding = "UTF-8")
}

sha256_file <- function(path) {
  if (!file.exists(path)) stop("Cannot hash missing file: ", path)
  output <- system2("sha256sum", path, stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) stop("sha256sum failed for ", path, ": ", paste(output, collapse = " "))
  hash <- strsplit(output[[1]], "[[:space:]]+")[[1]][[1]]
  if (!grepl("^[0-9a-f]{64}$", hash)) stop("Invalid SHA-256 output for ", path)
  hash
}

manifest_path <- file.path(discovery_root, "results", "05_report", "downstream_output_manifest.csv")
report_validation_path <- file.path(discovery_root, "results", "05_report", "05_report_validation.csv")
source_integrity_path <- file.path(discovery_root, "results", "05_report", "source_integrity_recheck.csv")
final_report_path <- file.path(discovery_root, "results", "05_report", "final_analysis_report.md")
audit_report_path <- file.path(audit_root, "cross_cohort_audit_report.md")
audit_validation_path <- file.path(audit_root, "00_cross_cohort_audit_validation.csv")
audit_pass_path <- file.path(cross_root, "state", "00_cross_cohort_audit.PASS")

manifest <- read_csv(manifest_path)
required_manifest_columns <- c("relative_path", "bytes", "sha256")
if (!identical(names(manifest), required_manifest_columns)) {
  stop("Discovery manifest schema changed: ", paste(names(manifest), collapse = ", "))
}
if (anyDuplicated(manifest$relative_path)) stop("Discovery manifest contains duplicate relative paths")
if (any(!grepl("^[0-9a-f]{64}$", manifest$sha256))) stop("Discovery manifest contains invalid SHA-256 values")
manifest$bytes <- as.numeric(manifest$bytes)
if (anyNA(manifest$bytes) || any(manifest$bytes < 0)) stop("Discovery manifest contains invalid byte counts")

report_validation <- read_csv(report_validation_path)
source_integrity <- read_csv(source_integrity_path)
audit_validation <- read_csv(audit_validation_path)
if (!file.exists(audit_pass_path)) stop("Cross-cohort audit PASS marker is missing")

manifest_sha256 <- sha256_file(manifest_path)
evidence_paths <- c(
  manifest_path,
  report_validation_path,
  source_integrity_path,
  final_report_path,
  audit_report_path,
  audit_validation_path,
  audit_pass_path
)
evidence_info <- file.info(evidence_paths)
freeze_evidence <- data.frame(
  evidence_type = c(
    "discovery_output_manifest",
    "discovery_final_validation",
    "discovery_source_integrity",
    "discovery_final_report",
    "cross_cohort_audit_report",
    "cross_cohort_audit_validation",
    "cross_cohort_audit_pass_marker"
  ),
  path = evidence_paths,
  bytes = as.numeric(evidence_info$size),
  sha256 = vapply(evidence_paths, sha256_file, character(1)),
  modified = format(evidence_info$mtime, "%Y-%m-%dT%H:%M:%S%z"),
  stringsAsFactors = FALSE
)
write_csv(freeze_evidence, "discovery_v1_freeze_evidence.csv")

key_relative_paths <- c(
  "results/03_annotate/10samples_annotated_seurat_v5.rds",
  "results/04_compare/all_cell_types_pseudobulk_de.csv.gz",
  "results/04_compare/all_cell_types_GO_BP_cameraPR.csv.gz",
  "results/05_report/final_analysis_report.md",
  "results/05_report/10samples_reproducibility_bundle.tar.gz"
)
key_assets <- manifest[match(key_relative_paths, manifest$relative_path), , drop = FALSE]
key_assets$present_in_manifest <- !is.na(key_assets$relative_path)
key_assets$relative_path[is.na(key_assets$relative_path)] <- key_relative_paths[is.na(key_assets$relative_path)]
write_csv(key_assets, "discovery_v1_key_assets.csv")

freeze <- data.frame(
  freeze_id = "discovery_10_v1",
  frozen_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  discovery_root = discovery_root,
  manifest_path = manifest_path,
  manifest_sha256 = manifest_sha256,
  manifest_entries = nrow(manifest),
  manifest_total_bytes = format(sum(manifest$bytes), scientific = FALSE, trim = TRUE),
  final_validation_pass = sum(report_validation$status == "PASS"),
  final_validation_total = nrow(report_validation),
  source_integrity_pass = sum(source_integrity$status == "PASS"),
  source_integrity_total = nrow(source_integrity),
  audit_validation_pass = sum(audit_validation$status == "PASS"),
  audit_validation_total = nrow(audit_validation),
  policy = "immutable_reference_only; do_not_overwrite; future_changes_require_new_version_root",
  stringsAsFactors = FALSE
)
write_csv(freeze, "discovery_v1_freeze_record.csv")

policy <- c(
  "# Discovery cohort v1 freeze",
  "",
  paste0("- Freeze ID: `", freeze$freeze_id, "`"),
  paste0("- Frozen at: ", freeze$frozen_at),
  paste0("- Source root: `", discovery_root, "`"),
  paste0("- Manifest entries: ", freeze$manifest_entries),
  paste0("- Manifest bytes: ", freeze$manifest_total_bytes),
  paste0("- Manifest SHA-256: `", manifest_sha256, "`"),
  paste0("- Final validation: ", freeze$final_validation_pass, "/", freeze$final_validation_total, " PASS"),
  paste0("- Source integrity: ", freeze$source_integrity_pass, "/", freeze$source_integrity_total, " PASS"),
  paste0("- Cross-cohort audit validation: ", freeze$audit_validation_pass, "/", freeze$audit_validation_total, " PASS"),
  "",
  "This freeze references the existing SHA-256 manifest. It does not copy or rehash every large result file.",
  "Validated discovery outputs must not be overwritten. Any later change requires a new versioned output root and a new manifest."
)
writeLines(policy, file.path(result, "discovery_v1_freeze_policy.md"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(result, "sessionInfo.txt"), useBytes = TRUE)

message("Discovery v1 freeze record created")
