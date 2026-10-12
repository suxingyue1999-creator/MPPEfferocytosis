#!/usr/bin/env Rscript

local_lib <- Sys.getenv("BD_R_LIB", "inputs/discovery/bd_rhapsody_10/downstream/10samples_scrna/env/Rlib")
.libPaths(c(local_lib, .libPaths()))
suppressPackageStartupMessages({library(Matrix); library(Seurat); library(ggplot2); library(digest)})
set.seed(20260726)
options(future.globals.maxSize = 40 * 1024^3)

root <- Sys.getenv("BD_ANALYSIS_ROOT", "inputs/discovery/bd_rhapsody_10/downstream/10samples_scrna")
input_root <- Sys.getenv("BD_INPUT_ROOT", "inputs/discovery/bd_rhapsody_10/work/outputs")
result <- file.path(root, "results", "03_annotate")
dir.create(result, recursive = TRUE, showWarnings = FALSE)
object <- readRDS(file.path(root, "results", "02_integrate", "10samples_integrated_seurat_v5.rds"))
DefaultAssay(object) <- "RNA"
object <- tryCatch(JoinLayers(object, assay = "RNA"), error = function(e) object)
object <- NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)

broad_signatures <- list(
  Macrophage_Monocyte = c("LST1", "TYROBP", "FCER1G", "CTSS", "AIF1", "LYZ", "CTSD", "LILRB1", "FCGR3A"),
  Neutrophil = c("FCGR3B", "CSF3R", "CEACAM8", "S100A8", "S100A9", "FPR1", "CXCR2", "MNDA"),
  T_NK = c("CD3D", "CD3E", "TRBC1", "TRBC2", "IL7R", "NKG7", "GNLY", "KLRD1"),
  B_Plasma = c("MS4A1", "CD79A", "CD79B", "CD37", "CD74", "MZB1", "JCHAIN", "SDC1"),
  Epithelial = c("EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "MUC1", "SFTPA1", "SFTPB"),
  Dendritic = c("FCER1A", "CD1C", "CLEC10A", "CD1E", "CLEC9A", "BATF3"),
  Endothelial = c("PECAM1", "VWF", "EMCN", "KDR", "ENG", "RAMP2"),
  Stromal = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "COL6A1"),
  Mast = c("KIT", "MS4A2", "TPSAB1", "TPSB2", "CPA3", "HDC"),
  Platelet = c("PPBP", "PF4", "NRGN", "GP9", "ITGA2B")
)

restore_average_group_names <- function(observed, expected) {
  observed <- as.character(observed)
  expected <- unique(as.character(expected))
  if (!length(expected) || anyNA(expected) || any(!nzchar(expected))) {
    stop("Grouping metadata contains missing or empty values")
  }
  canonical_expected <- gsub("_", "-", expected, fixed = TRUE)
  numeric_start <- grepl("^[0-9]", canonical_expected)
  canonical_expected[numeric_start] <- paste0("g", canonical_expected[numeric_start])
  mapped <- vapply(observed, function(value) {
    exact <- which(expected == value)
    if (length(exact) == 1L) return(expected[[exact]])
    canonical <- which(canonical_expected == value)
    if (length(canonical) == 1L) return(expected[[canonical]])
    NA_character_
  }, character(1), USE.NAMES = FALSE)
  if (any(!mapped %in% expected) || anyDuplicated(mapped) || !setequal(mapped, expected)) {
    stop(
      "Could not map AverageExpression groups to metadata values: observed=",
      paste(observed, collapse = ","), "; expected=", paste(expected, collapse = ",")
    )
  }
  mapped
}

score_groups <- function(seurat_object, group_column, signatures) {
  group_values <- seurat_object[[]][[group_column]]
  if (is.null(group_values)) stop("Missing grouping metadata column: ", group_column)
  average <- AverageExpression(seurat_object, assays = "RNA", layer = "data", group.by = group_column, verbose = FALSE)$RNA
  groups <- restore_average_group_names(colnames(average), group_values)
  colnames(average) <- groups
  z <- t(scale(t(as.matrix(average))))
  z[!is.finite(z)] <- 0
  rows <- list()
  for (label in names(signatures)) {
    genes <- intersect(signatures[[label]], rownames(average))
    score <- if (length(genes)) colMeans(z[genes, , drop = FALSE]) else rep(-Inf, length(groups))
    support <- if (length(genes)) colSums(average[genes, , drop = FALSE] > 0.1) else rep(0, length(groups))
    rows[[label]] <- data.frame(group = groups, candidate = label, score = score, supporting_markers = support)
  }
  scores <- do.call(rbind, rows)
  assignments <- lapply(split(scores, scores$group), function(frame) {
    frame <- frame[order(frame$score, decreasing = TRUE), , drop = FALSE]
    margin <- frame$score[1] - frame$score[2]
    label <- if (frame$score[1] > 0 && frame$supporting_markers[1] >= 2) frame$candidate[1] else "Unresolved"
    data.frame(
      group = frame$group[1], label = label, top_score = frame$score[1],
      score_margin = margin, supporting_markers = frame$supporting_markers[1],
      confidence = ifelse(label == "Unresolved", "low", ifelse(margin >= 0.35, "high", "moderate"))
    )
  })
  list(scores = scores, assignments = do.call(rbind, assignments))
}

message("Finding broad cluster markers")
Idents(object) <- object$seurat_clusters
marker_workers <- suppressWarnings(as.integer(Sys.getenv("BD_MARKER_WORKERS", "12")))
if (!is.finite(marker_workers) || marker_workers < 1L) marker_workers <- 1L
if (marker_workers > 1L) future::plan(future::multicore, workers = marker_workers)
markers <- FindAllMarkers(
  object, assay = "RNA", only.pos = TRUE, test.use = "wilcox",
  min.pct = 0.20, logfc.threshold = 0.25, max.cells.per.ident = 3000,
  random.seed = 20260726, verbose = FALSE
)
future::plan(future::sequential)
write.csv(markers, gzfile(file.path(result, "broad_cluster_markers.csv.gz")), row.names = FALSE)

broad <- score_groups(object, "seurat_clusters", broad_signatures)
colnames(broad$assignments)[colnames(broad$assignments) == "group"] <- "cluster"
colnames(broad$scores)[colnames(broad$scores) == "group"] <- "cluster"
write.csv(broad$scores, file.path(result, "broad_annotation_scores.csv"), row.names = FALSE)
write.csv(broad$assignments, file.path(result, "broad_cluster_annotations.csv"), row.names = FALSE)
map <- setNames(broad$assignments$label, broad$assignments$cluster)
confidence_map <- setNames(broad$assignments$confidence, broad$assignments$cluster)
object$cell_type_broad <- unname(map[as.character(object$seurat_clusters)])
object$annotation_confidence <- unname(confidence_map[as.character(object$seurat_clusters)])
if (anyNA(object$cell_type_broad) || any(!nzchar(object$cell_type_broad))) {
  stop("Broad annotation mapping produced missing labels")
}
if (anyNA(object$annotation_confidence) || any(!nzchar(object$annotation_confidence))) {
  stop("Broad annotation mapping produced missing confidence values")
}
object$cell_type_subtype <- object$cell_type_broad

subtype_signatures <- list(
  Macrophage_Monocyte = list(
    Classical_Monocyte = c("S100A8", "S100A9", "FCN1", "CTSD", "LYZ", "VCAN"),
    Nonclassical_Monocyte = c("FCGR3A", "MS4A7", "LST1", "LILRB1", "IFITM3"),
    Alveolar_Macrophage = c("FABP4", "PPARG", "MARCO", "INHBA", "MSR1"),
    Inflammatory_Macrophage = c("SPP1", "APOE", "CTSB", "CCL3", "CCL4", "IL1B"),
    Cycling_Myeloid = c("MKI67", "TOP2A", "TYMS", "STMN1")
  ),
  Neutrophil = list(
    Mature_Neutrophil = c("FCGR3B", "CSF3R", "CXCR2", "FPR1", "CEACAM8"),
    Activated_Neutrophil = c("S100A8", "S100A9", "FOS", "JUN", "IL1B"),
    IFN_Neutrophil = c("ISG15", "IFIT1", "IFIT3", "MX1", "IFI6")
  ),
  T_NK = list(
    CD4_T = c("CD3D", "CD4", "IL7R", "LTB", "MAL"),
    CD8_T = c("CD3D", "CD8A", "CD8B", "CCL5", "GZMK"),
    NK = c("NKG7", "GNLY", "KLRD1", "PRF1", "GZMB"),
    Regulatory_T = c("FOXP3", "IL2RA", "CTLA4", "TIGIT"),
    Cycling_T_NK = c("MKI67", "TOP2A", "TYMS", "STMN1")
  ),
  B_Plasma = list(
    Naive_B = c("MS4A1", "TCL1A", "IGHD", "IL4R", "FCER2"),
    Memory_B = c("CD27", "CD37", "AIM2", "BANK1", "CD79A"),
    Plasma = c("MZB1", "JCHAIN", "SDC1", "XBP1", "IGHG1")
  ),
  Epithelial = list(
    AT2 = c("SFTPA1", "SFTPA2", "SFTPB", "SFTPC", "ABCA3"),
    AT1 = c("AGER", "CAV1", "HOPX", "CAV2", "EMP2"),
    Ciliated = c("FOXJ1", "PIFO", "TPPP3", "CAPS", "TUBA1A"),
    Secretory = c("SCGB1A1", "SCGB3A1", "MUC1", "KRT19"),
    Basal = c("KRT5", "KRT14", "TP63", "KRT17"),
    Cycling_Epithelial = c("MKI67", "TOP2A", "TYMS", "STMN1")
  )
)

lineage_workers <- suppressWarnings(as.integer(Sys.getenv("BD_LINEAGE_WORKERS", "5")))
if (!is.finite(lineage_workers) || lineage_workers < 1L) lineage_workers <- 1L
lineage_workers <- min(lineage_workers, length(subtype_signatures))

recluster_lineage <- function(lineage) {
  cell_index <- which(!is.na(object$cell_type_broad) & object$cell_type_broad == lineage)
  cells <- colnames(object)[cell_index]
  if (length(cells) < 500) {
    message("Skipping lineage with <500 cells: ", lineage)
    return(NULL)
  }
  set.seed(20260726 + match(lineage, names(subtype_signatures)))
  message("Reclustering lineage: ", lineage, " (", length(cells), " cells)")
  subset_object <- subset(object, cells = cells)
  DefaultAssay(subset_object) <- "RNA"
  subset_object <- FindVariableFeatures(subset_object, nfeatures = 2000, verbose = FALSE)
  subset_object <- ScaleData(subset_object, features = VariableFeatures(subset_object), verbose = FALSE)
  subset_object <- RunPCA(subset_object, features = VariableFeatures(subset_object), npcs = 30, approx = TRUE, verbose = FALSE)
  subset_object <- FindNeighbors(subset_object, dims = 1:20, verbose = FALSE)
  subset_object <- FindClusters(subset_object, resolution = 0.4, algorithm = 1, random.seed = 20260726, verbose = FALSE)
  subset_object$lineage_cluster <- paste(lineage, as.character(Idents(subset_object)), sep = "_")
  subset_object <- RunUMAP(subset_object, dims = 1:20, seed.use = 20260726, verbose = FALSE)
  scored <- score_groups(subset_object, "lineage_cluster", subtype_signatures[[lineage]])
  subtype_map <- setNames(scored$assignments$label, scored$assignments$group)
  subtype <- unname(subtype_map[subset_object$lineage_cluster])
  if (anyNA(subtype) || any(!nzchar(subtype))) stop("Subtype mapping produced missing labels for lineage: ", lineage)
  subtype[subtype == "Unresolved"] <- paste0(lineage, "_Unresolved")
  emb <- as.data.frame(Embeddings(subset_object, "umap"))
  emb$cell_id <- rownames(emb)
  emb$lineage <- lineage
  emb$lineage_cluster <- subset_object$lineage_cluster
  emb$cell_type_subtype <- subtype
  assignment <- scored$assignments
  assignment$lineage <- lineage
  p <- DimPlot(subset_object, reduction = "umap", group.by = "lineage_cluster", label = TRUE, repel = TRUE, raster = TRUE) + ggtitle(lineage) + NoLegend()
  ggsave(file.path(result, paste0("lineage_", lineage, "_umap.pdf")), p, width = 8, height = 6)
  output <- list(
    lineage = lineage,
    cell_ids = colnames(subset_object),
    subtype = subtype,
    embedding = emb,
    assignment = assignment
  )
  rm(subset_object, p)
  gc(verbose = FALSE)
  output
}

message("Lineage reclustering workers: ", lineage_workers)
lineage_results <- parallel::mclapply(
  names(subtype_signatures),
  recluster_lineage,
  mc.cores = lineage_workers,
  mc.preschedule = FALSE,
  mc.set.seed = TRUE
)
failed_lineages <- vapply(lineage_results, inherits, logical(1), what = "try-error")
if (any(failed_lineages)) {
  failure_report <- data.frame(
    lineage = names(subtype_signatures)[failed_lineages],
    error = vapply(lineage_results[failed_lineages], as.character, character(1)),
    stringsAsFactors = FALSE
  )
  write.csv(failure_report, file.path(result, "lineage_worker_failures.csv"), row.names = FALSE)
  stop("Lineage worker failure: ", paste(failure_report$lineage, collapse = ","))
}
lineage_results <- Filter(Negate(is.null), lineage_results)
lineage_embeddings <- list()
lineage_assignments <- list()
for (lineage_result in lineage_results) {
  lineage <- lineage_result$lineage
  object$cell_type_subtype[match(lineage_result$cell_ids, colnames(object))] <- lineage_result$subtype
  lineage_embeddings[[lineage]] <- lineage_result$embedding
  lineage_assignments[[lineage]] <- lineage_result$assignment
}
write.csv(do.call(rbind, lineage_embeddings), gzfile(file.path(result, "lineage_recluster_embeddings.csv.gz")), row.names = FALSE)
write.csv(do.call(rbind, lineage_assignments), file.path(result, "lineage_subcluster_annotations.csv"), row.names = FALSE)

bd_rows <- list()
label_manifest <- list()
for (sample_id in unique(object$sample_id)) {
  path <- file.path(input_root, sample_id, paste0(sample_id, "_cell_type_experimental.csv"))
  labels <- read.csv(path, stringsAsFactors = FALSE)
  bd_rows[[sample_id]] <- data.frame(
    cell_id = paste(sample_id, labels$Cell_Index, sep = "_"),
    bd_experimental_label = labels$Cell_Type_Experimental,
    stringsAsFactors = FALSE
  )
  info <- file.info(path)
  label_manifest[[sample_id]] <- data.frame(sample_id = sample_id, path = path, bytes = info$size, sha256 = digest(path, algo = "sha256", file = TRUE))
}
bd <- do.call(rbind, bd_rows)
label_index <- match(colnames(object), bd$cell_id)
object$bd_experimental_label <- bd$bd_experimental_label[label_index]
write.csv(do.call(rbind, label_manifest), file.path(result, "bd_experimental_label_manifest.csv"), row.names = FALSE)

map_bd <- function(values) {
  output <- rep("Other", length(values))
  output[grepl("^(T_|T$|NK)", values, ignore.case = TRUE)] <- "T_NK"
  output[grepl("^(B$|B_|Plasma)", values, ignore.case = TRUE)] <- "B_Plasma"
  output[grepl("Neut", values, ignore.case = TRUE)] <- "Neutrophil"
  output[grepl("Mono|Macro", values, ignore.case = TRUE)] <- "Macrophage_Monocyte"
  output[grepl("Dend", values, ignore.case = TRUE)] <- "Dendritic"
  output[grepl("Epit|AT1|AT2|Cilia|Secret|Basal", values, ignore.case = TRUE)] <- "Epithelial"
  output[grepl("Endothel", values, ignore.case = TRUE)] <- "Endothelial"
  output[grepl("Fibro|Strom", values, ignore.case = TRUE)] <- "Stromal"
  output[grepl("Mast", values, ignore.case = TRUE)] <- "Mast"
  output[grepl("Platelet", values, ignore.case = TRUE)] <- "Platelet"
  output
}
object$bd_experimental_broad <- map_bd(object$bd_experimental_label)
concordance <- as.data.frame(table(final_broad = object$cell_type_broad, bd_reference_broad = object$bd_experimental_broad), stringsAsFactors = FALSE)
concordance <- concordance[concordance$Freq > 0, , drop = FALSE]
write.csv(concordance, file.path(result, "bd_reference_concordance.csv"), row.names = FALSE)

annotations <- object[[]][, c("sample_id", "severity", "barcode", "seurat_clusters", "cell_type_broad", "cell_type_subtype", "annotation_confidence", "bd_experimental_label")]
annotations$cell_id <- rownames(annotations)
annotations <- annotations[, c("cell_id", setdiff(colnames(annotations), "cell_id"))]
write.csv(annotations, gzfile(file.path(result, "cell_annotations.csv.gz")), row.names = FALSE)
write.csv(as.data.frame(table(cell_type_broad = object$cell_type_broad, sample_id = object$sample_id), stringsAsFactors = FALSE), file.path(result, "cell_type_counts_by_sample.csv"), row.names = FALSE)
saveRDS(object, file.path(result, "10samples_annotated_seurat_v5.rds"), compress = FALSE)

p_broad <- DimPlot(object, reduction = "umap", group.by = "cell_type_broad", label = TRUE, repel = TRUE, raster = TRUE)
p_subtype <- DimPlot(object, reduction = "umap", group.by = "cell_type_subtype", label = FALSE, raster = TRUE)
ggsave(file.path(result, "umap_broad_cell_types.pdf"), p_broad, width = 11, height = 8)
ggsave(file.path(result, "umap_broad_cell_types.png"), p_broad, width = 11, height = 8, dpi = 170)
ggsave(file.path(result, "umap_cell_subtypes.pdf"), p_subtype, width = 12, height = 8)

writeLines(c(
  "# Stage 03: cell annotation",
  "",
  "- Final broad labels were assigned from cluster-level expression of curated canonical marker signatures.",
  paste0("- Parallelism: ", marker_workers, " marker workers and ", lineage_workers, " independent lineage workers."),
  "- Major lineages with at least 500 cells were independently normalized, reclustered, and assigned subtype labels.",
  "- Low-evidence clusters remain explicitly unresolved; no forced label is used.",
  "- BD experimental labels were SHA-256 inventoried and used only for a post hoc concordance table.",
  "- BD experimental labels did not enter marker scoring, clustering, or final label assignment.",
  "- Final biological interpretation remains subject to marker-plot review."
), file.path(result, "annotation_report.md"))
writeLines(capture.output(sessionInfo()), file.path(result, "r_sessionInfo.txt"))
