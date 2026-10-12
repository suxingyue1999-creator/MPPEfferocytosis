#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

cross_root <- Sys.getenv(
  "BD_CROSS_ROOT",
  "inputs/discovery/bd_rhapsody_10/downstream/10_23_cross_cohort_v1"
)
result <- file.path(cross_root, "results", "01_discovery_freeze")

read_csv <- function(name) {
  path <- file.path(result, name)
  if (!file.exists(path)) stop("Missing freeze output: ", path)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8-BOM")
}

sha256_file <- function(path) {
  output <- system2("sha256sum", path, stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) return(NA_character_)
  strsplit(output[[1]], "[[:space:]]+")[[1]][[1]]
}

freeze <- read_csv("discovery_v1_freeze_record.csv")
evidence <- read_csv("discovery_v1_freeze_evidence.csv")
key_assets <- read_csv("discovery_v1_key_assets.csv")
policy_path <- file.path(result, "discovery_v1_freeze_policy.md")

checks <- list()
add_check <- function(check, status, detail) {
  checks[[length(checks) + 1L]] <<- data.frame(
    check = check,
    status = if (isTRUE(status)) "PASS" else "FAIL",
    detail = as.character(detail),
    stringsAsFactors = FALSE
  )
}

add_check("single_freeze_record", nrow(freeze) == 1L, paste0("rows=", nrow(freeze)))
add_check("freeze_id_locked", identical(freeze$freeze_id, "discovery_10_v1"), freeze$freeze_id)
add_check(
  "manifest_completion_totals",
  as.integer(freeze$manifest_entries) == 215L && as.numeric(freeze$manifest_total_bytes) == 33038526565,
  paste0("entries=", freeze$manifest_entries, ";bytes=", freeze$manifest_total_bytes)
)
add_check(
  "final_validation_7_of_7",
  as.integer(freeze$final_validation_pass) == 7L && as.integer(freeze$final_validation_total) == 7L,
  paste0(freeze$final_validation_pass, "/", freeze$final_validation_total)
)
add_check(
  "source_integrity_60_of_60",
  as.integer(freeze$source_integrity_pass) == 60L && as.integer(freeze$source_integrity_total) == 60L,
  paste0(freeze$source_integrity_pass, "/", freeze$source_integrity_total)
)
add_check(
  "audit_validation_all_pass",
  as.integer(freeze$audit_validation_pass) == as.integer(freeze$audit_validation_total) && as.integer(freeze$audit_validation_total) > 0L,
  paste0(freeze$audit_validation_pass, "/", freeze$audit_validation_total)
)

evidence_exists <- file.exists(evidence$path)
current_hash <- rep(NA_character_, nrow(evidence))
current_hash[evidence_exists] <- vapply(evidence$path[evidence_exists], sha256_file, character(1))
hash_match <- evidence_exists & !is.na(current_hash) & current_hash == evidence$sha256
add_check(
  "freeze_evidence_hashes_match",
  all(hash_match),
  paste0("files=", nrow(evidence), ";matching=", sum(hash_match))
)
manifest_row <- evidence$evidence_type == "discovery_output_manifest"
add_check(
  "manifest_hash_locked",
  sum(manifest_row) == 1L && identical(freeze$manifest_sha256, evidence$sha256[manifest_row]),
  freeze$manifest_sha256
)
add_check(
  "key_assets_present_in_manifest",
  nrow(key_assets) == 5L && all(key_assets$present_in_manifest %in% c(TRUE, "TRUE")) && all(grepl("^[0-9a-f]{64}$", key_assets$sha256)),
  paste0("assets=", nrow(key_assets))
)
add_check(
  "immutable_policy_present",
  file.exists(policy_path) && grepl("must not be overwritten", paste(readLines(policy_path, warn = FALSE), collapse = " "), fixed = TRUE),
  policy_path
)

validation <- do.call(rbind, checks)
write.csv(validation, file.path(result, "01_discovery_freeze_validation.csv"), row.names = FALSE)
report <- c(
  "# Discovery v1 freeze validation",
  "",
  capture.output(print(validation, row.names = FALSE))
)
writeLines(report, file.path(result, "01_discovery_freeze_validation.md"), useBytes = TRUE)

if (any(validation$status != "PASS")) {
  stop("Discovery v1 freeze validation failed: ", paste(validation$check[validation$status != "PASS"], collapse = ", "))
}
message("Discovery v1 freeze validation PASS")
