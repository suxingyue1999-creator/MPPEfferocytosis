#!/usr/bin/env Rscript

lib <- Sys.getenv("BD_R_LIB", "inputs/discovery/bd_rhapsody_10/downstream/10samples_scrna/env/Rlib")
.libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages(library(Seurat))
root <- Sys.getenv("BD_ANALYSIS_ROOT", "inputs/discovery/bd_rhapsody_10/downstream/10samples_scrna")
result <- file.path(root, "results", "03_annotate")
object <- readRDS(file.path(result, "10samples_annotated_seurat_v5.rds"))
checks <- list()
add <- function(check, passed, detail) checks[[length(checks) + 1L]] <<- data.frame(check = check, status = if (passed) "PASS" else "FAIL", detail = detail)
add("all_cells_have_broad_label", all(!is.na(object$cell_type_broad) & nzchar(object$cell_type_broad)), paste0("cells=", ncol(object)))
add("all_cells_have_subtype_label", all(!is.na(object$cell_type_subtype) & nzchar(object$cell_type_subtype)), paste0("cells=", ncol(object)))
core <- c("Macrophage_Monocyte", "Neutrophil", "T_NK", "B_Plasma", "Epithelial")
allowed_broad <- c(core, "Dendritic", "Endothelial", "Stromal", "Mast", "Platelet", "Unresolved")
unexpected_broad <- setdiff(unique(object$cell_type_broad), allowed_broad)
add("broad_label_domain", !length(unexpected_broad), paste(unexpected_broad, collapse = ","))
present <- intersect(core, unique(object$cell_type_broad))
add("core_lineages_detected", length(present) >= 4, paste(present, collapse = ","))
cluster_labels <- aggregate(cell_type_broad ~ seurat_clusters, object[[]], function(x) length(unique(x)))
add("broad_labels_cluster_consistent", all(cluster_labels$cell_type_broad == 1L), paste0("clusters=", nrow(cluster_labels)))
broad_clusters <- read.csv(file.path(result, "broad_cluster_annotations.csv"), stringsAsFactors = FALSE)$cluster
object_clusters <- unique(as.character(object$seurat_clusters))
add(
  "broad_annotation_clusters_match",
  !anyNA(broad_clusters) && !anyDuplicated(broad_clusters) && setequal(broad_clusters, object_clusters),
  paste0("annotation=", length(broad_clusters), "; object=", length(object_clusters))
)
reclustered <- read.csv(file.path(result, "lineage_subcluster_annotations.csv"), stringsAsFactors = FALSE)
add("major_lineage_reclustering", length(unique(reclustered$lineage)) >= 3, paste(unique(reclustered$lineage), collapse = ","))
lineage_embeddings <- read.csv(gzfile(file.path(result, "lineage_recluster_embeddings.csv.gz")), stringsAsFactors = FALSE)
broad_counts <- table(object$cell_type_broad)
expected_reclustered <- intersect(core, names(broad_counts)[broad_counts >= 500])
actual_reclustered <- unique(lineage_embeddings$lineage)
add(
  "eligible_lineages_reclustered",
  !anyDuplicated(lineage_embeddings$cell_id) && setequal(actual_reclustered, expected_reclustered),
  paste0("expected=", paste(expected_reclustered, collapse = ","), "; actual=", paste(actual_reclustered, collapse = ","))
)
embedding_index <- match(lineage_embeddings$cell_id, colnames(object))
embedding_labels_match <- !anyNA(embedding_index) && all(object$cell_type_broad[embedding_index] == lineage_embeddings$lineage)
add("lineage_cells_match_broad_labels", embedding_labels_match, paste0("cells=", nrow(lineage_embeddings)))
manifest <- read.csv(file.path(result, "bd_experimental_label_manifest.csv"), stringsAsFactors = FALSE)
add("bd_reference_manifest_10", nrow(manifest) == 10 && all(nchar(manifest$sha256) == 64), paste0("rows=", nrow(manifest)))
required <- c("broad_cluster_markers.csv.gz", "broad_annotation_scores.csv", "broad_cluster_annotations.csv", "lineage_recluster_embeddings.csv.gz", "lineage_subcluster_annotations.csv", "cell_annotations.csv.gz", "bd_reference_concordance.csv", "umap_broad_cell_types.pdf", "annotation_report.md")
missing <- required[!file.exists(file.path(result, required)) | file.info(file.path(result, required))$size == 0]
add("required_outputs", !length(missing), paste(missing, collapse = ","))
add("sample_metadata_preserved", length(unique(object$sample_id)) == 10 && length(unique(object$severity)) == 2, paste(unique(object$severity), collapse = ","))
validation <- do.call(rbind, checks)
write.csv(validation, file.path(result, "03_annotate_validation.csv"), row.names = FALSE)
writeLines(c("# Stage 03 validation", "", paste0("- ", validation$status, ": ", validation$check, " (", validation$detail, ")")), file.path(result, "03_annotate_validation.md"))
print(validation, row.names = FALSE)
if (any(validation$status != "PASS")) quit(status = 1L)
