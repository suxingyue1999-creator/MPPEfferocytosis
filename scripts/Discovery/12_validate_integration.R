#!/usr/bin/env Rscript

local_lib <- Sys.getenv("BD_R_LIB", "inputs/discovery/bd_rhapsody_10/downstream/10samples_scrna/env/Rlib")
.libPaths(c(local_lib, .libPaths()))
suppressPackageStartupMessages({library(Seurat); library(Matrix)})
root <- Sys.getenv("BD_ANALYSIS_ROOT", "inputs/discovery/bd_rhapsody_10/downstream/10samples_scrna")
result <- file.path(root, "results", "02_integrate")
object <- readRDS(file.path(result, "10samples_integrated_seurat_v5.rds"))
keep <- read.csv(gzfile(file.path(root, "results", "01_qc", "cells_keep.csv.gz")), stringsAsFactors = FALSE)
samples <- c("D01", "D02", "D03", "D04", "D05", "D06", "D07", "D08", "D09", "D10")

checks <- list()
add <- function(check, passed, detail) {
  checks[[length(checks) + 1L]] <<- data.frame(check = check, status = if (passed) "PASS" else "FAIL", detail = detail)
}
add("seurat_v5", packageVersion("Seurat") >= "5.0.0", paste0("version=", packageVersion("Seurat")))
add("cell_count_matches_qc", ncol(object) == nrow(keep), paste0("object=", ncol(object), "; expected=", nrow(keep)))
add("cell_ids_unique", !anyDuplicated(colnames(object)), paste0("unique=", length(unique(colnames(object)))))
add("sample_ids_fixed", identical(sort(unique(object$sample_id)), sort(samples)), paste(unique(object$sample_id), collapse = ","))
sample_severity <- unique(object[[]][, c("sample_id", "severity")])
add("severity_5v5", identical(unname(as.integer(sort(table(sample_severity$severity)))), c(5L, 5L)), paste(capture.output(print(table(sample_severity$severity))), collapse = " "))
add("required_assays", all(c("RNA", "SCT", "integrated") %in% Assays(object)), paste(Assays(object), collapse = ","))
add("required_reductions", all(c("pca", "umap") %in% Reductions(object)), paste(Reductions(object), collapse = ","))
clusters <- unique(object$seurat_clusters)
add("cluster_guardrails", length(clusters) >= 5 && length(clusters) <= 100, paste0("clusters=", length(clusters)))
umap <- Embeddings(object, "umap")
add("umap_finite", nrow(umap) == ncol(object) && all(is.finite(umap)), paste(dim(umap), collapse = "x"))
DefaultAssay(object) <- "RNA"
counts <- tryCatch(LayerData(object, assay = "RNA", layer = "counts"), error = function(e) NULL)
if (is.null(counts)) counts <- GetAssayData(object, assay = "RNA", slot = "counts")
probe <- counts@x[seq_len(min(length(counts@x), 1000000L))]
add("raw_counts_present_integer", ncol(counts) == ncol(object) && all(probe >= 0) && all(probe == round(probe)), paste(dim(counts), collapse = "x"))

validation <- do.call(rbind, checks)
write.csv(validation, file.path(result, "02_integrate_validation.csv"), row.names = FALSE)
writeLines(c("# Stage 02 validation", "", paste0("- ", validation$status, ": ", validation$check, " (", validation$detail, ")")), file.path(result, "02_integrate_validation.md"))
print(validation, row.names = FALSE)
if (any(validation$status != "PASS")) quit(status = 1L)
