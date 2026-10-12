#!/usr/bin/env Rscript

cross_root <- Sys.getenv(
  "BD_CROSS_ROOT",
  "inputs/discovery/bd_rhapsody_10/downstream/10_23_cross_cohort_v1"
)
result <- file.path(cross_root, "results", "00_audit")

read_csv <- function(name) {
  path <- file.path(result, name)
  if (!file.exists(path)) stop("Missing audit output: ", path)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8-BOM")
}

checks <- list()
add <- function(check, status, detail = "") {
  checks[[length(checks) + 1L]] <<- data.frame(
    check = check,
    status = if (isTRUE(status)) "PASS" else "FAIL",
    detail = detail,
    stringsAsFactors = FALSE
  )
}

inventory <- read_csv("cohort_inventory.csv")
add(
  "cohort_inventory_exact",
  identical(inventory$cohort, c("discovery_10", "validation_23")) &&
    identical(as.integer(inventory$samples), c(10L, 23L)) &&
    identical(as.integer(inventory$cells), c(115594L, 209299L)),
  paste(inventory$cohort, inventory$samples, inventory$cells, collapse = ";")
)

source_checks <- read_csv("source_validation_checks.csv")
add(
  "discovery_source_validators_pass",
  nrow(source_checks) == 27L && all(source_checks$status == "PASS"),
  paste0("rows=", nrow(source_checks), ";failures=", sum(source_checks$status != "PASS"))
)

broad <- read_csv("discovery_broad_cluster_review.csv")
add(
  "discovery_cluster_review_complete",
  nrow(broad) == 27L && all(broad$classification %in% c("通过", "需修订", "暂不可信")),
  paste0("clusters=", nrow(broad))
)

composition <- read_csv("discovery_composition_loo_review.csv")
add(
  "composition_leave_one_out_complete",
  nrow(composition) == 10L && all(is.finite(composition$full_severe_minus_mild)) &&
    all(is.finite(composition$dominant_cell_share)),
  paste0("cell_types=", nrow(composition))
)

de_review <- read_csv("discovery_pseudobulk_review.csv")
add(
  "discovery_pseudobulk_review_complete",
  nrow(de_review) == 7L && all(de_review$classification %in% c("通过", "需修订", "暂不可信")) &&
    de_review$classification[de_review$cell_type_broad == "Unresolved"] == "暂不可信",
  paste0("cell_types=", nrow(de_review))
)

validation <- read_csv("validation_evidence_review.csv")
add(
  "validation_evidence_review_complete",
  nrow(validation) == 9L && all(validation$classification %in% c("通过", "需修订", "暂不可信")),
  paste0("items=", nrow(validation))
)

mapping <- read_csv("cross_cohort_celltype_mapping.csv")
add(
  "harmonized_mapping_prespecified",
  sum(mapping$meta_eligible) == 4L &&
    all(c("B_lineage", "Macrophage_Monocyte", "Neutrophil", "T_lineage") %in% mapping$harmonized_label),
  paste(mapping$harmonized_label[mapping$meta_eligible], collapse = ",")
)

inputs <- read_csv("input_inventory.csv")
add(
  "required_inputs_present",
  nrow(inputs) == 6L && all(inputs$exists) && all(inputs$bytes > 0),
  paste0("files=", nrow(inputs), ";missing=", sum(!inputs$exists))
)

required <- c(
  "classification_summary.csv",
  "cross_cohort_audit_report.md",
  "discovery_lineage_subcluster_review.csv",
  "sessionInfo.txt"
)
add(
  "required_audit_outputs",
  all(file.exists(file.path(result, required))),
  paste(required[!file.exists(file.path(result, required))], collapse = ",")
)

validation_table <- do.call(rbind, checks)
write.csv(
  validation_table,
  file.path(result, "00_cross_cohort_audit_validation.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
writeLines(
  c(
    "# Cross-cohort audit validation",
    "",
    capture.output(print(validation_table, row.names = FALSE))
  ),
  file.path(result, "00_cross_cohort_audit_validation.md"),
  useBytes = TRUE
)

if (any(validation_table$status != "PASS")) {
  stop("Cross-cohort audit validation failed: ", paste(validation_table$check[validation_table$status != "PASS"], collapse = ","))
}
message("Cross-cohort audit validation PASS")
