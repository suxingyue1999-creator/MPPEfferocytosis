#!/usr/bin/env Rscript

local_lib <- Sys.getenv("BD_R_LIB", "inputs/discovery/bd_rhapsody_10/downstream/10samples_scrna/env/Rlib")
.libPaths(c(local_lib, .libPaths()))
suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(ggplot2)
})

set.seed(20260726)
options(future.globals.maxSize = 40 * 1024^3)
samples <- c("D01", "D02", "D03", "D04", "D05", "D06", "D07", "D08", "D09", "D10")
input_root <- Sys.getenv("BD_INPUT_ROOT", "inputs/discovery/bd_rhapsody_10/work/outputs")
analysis_root <- Sys.getenv("BD_ANALYSIS_ROOT", "inputs/discovery/bd_rhapsody_10/downstream/10samples_scrna")
result_dir <- file.path(analysis_root, "results", "02_integrate")
tmp_root <- file.path(analysis_root, "tmp", "02_integrate")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tmp_root, recursive = TRUE, showWarnings = FALSE)

locked <- read.csv(file.path(analysis_root, "results", "00_audit", "10samples_metadata_locked.csv"), stringsAsFactors = FALSE)
qc <- read.csv(gzfile(file.path(analysis_root, "results", "01_qc", "cell_qc.csv.gz")), stringsAsFactors = FALSE)
qc <- qc[qc$qc_keep, , drop = FALSE]

read_filtered_mex <- function(sample_id) {
  archive <- file.path(input_root, sample_id, paste0(sample_id, "_RSEC_MolsPerCell_MEX.zip"))
  target <- file.path(tmp_root, sample_id)
  unlink(target, recursive = TRUE, force = TRUE)
  dir.create(target, recursive = TRUE, showWarnings = FALSE)
  utils::unzip(archive, exdir = target)
  value <- Read10X(target, gene.column = 2, cell.column = 1, unique.features = TRUE)
  if (is.list(value)) value <- value[["Gene Expression"]]
  as(value, "dgCMatrix")
}

sample_workers <- suppressWarnings(as.integer(Sys.getenv("BD_SCT_WORKERS", "3")))
if (!is.finite(sample_workers) || sample_workers < 1L) sample_workers <- 1L
sample_workers <- min(sample_workers, length(samples))

transform_sample <- function(index) {
  sample_id <- samples[[index]]
  set.seed(20260726 + index)
  message("SCTransform: ", sample_id)
  counts <- read_filtered_mex(sample_id)
  sample_qc <- qc[qc$sample_id == sample_id, , drop = FALSE]
  index <- match(sample_qc$barcode, colnames(counts))
  if (anyNA(index)) stop(sample_id, ": retained cell absent from canonical filtered MEX")
  counts <- counts[, index, drop = FALSE]
  colnames(counts) <- sample_qc$cell_id
  object <- CreateSeuratObject(counts = counts, project = sample_id, min.cells = 0, min.features = 0)
  rownames(sample_qc) <- sample_qc$cell_id
  object <- AddMetaData(object, sample_qc[colnames(object), c(
    "sample_id", "severity", "barcode", "nCount_RNA", "nFeature_RNA",
    "percent_mt", "percent_ribo", "percent_hb", "doublet_score",
    "ambient_top20_percent"
  )])
  object$cohort_id <- locked$cohort_id[match(sample_id, locked$sample_id)]
  object <- SCTransform(
    object,
    assay = "RNA",
    new.assay.name = "SCT",
    method = "glmGamPoi",
    vst.flavor = "v2",
    vars.to.regress = "percent_mt",
    variable.features.n = 3000,
    conserve.memory = TRUE,
    return.only.var.genes = TRUE,
    seed.use = 20260726,
    verbose = TRUE
  )
  rm(counts)
  unlink(file.path(tmp_root, sample_id), recursive = TRUE, force = TRUE)
  gc(verbose = FALSE)
  object
}

message("SCTransform sample workers: ", sample_workers)
objects <- parallel::mclapply(
  seq_along(samples),
  transform_sample,
  mc.cores = sample_workers,
  mc.preschedule = FALSE,
  mc.set.seed = TRUE
)
failed_workers <- vapply(objects, inherits, logical(1), what = "try-error")
if (any(failed_workers)) stop("SCTransform sample worker failure: ", paste(samples[failed_workers], collapse = ","))
names(objects) <- samples

message("Selecting SCT integration features")
features <- SelectIntegrationFeatures(object.list = objects, nfeatures = 3000)
objects <- PrepSCTIntegration(object.list = objects, anchor.features = features, verbose = TRUE)
objects <- lapply(objects, function(object) {
  RunPCA(object, features = features, npcs = 40, approx = TRUE, verbose = FALSE)
})
cell_counts <- vapply(objects, ncol, integer(1))
reference_index <- which.min(abs(cell_counts - median(cell_counts)))
message("RPCA reference: ", names(objects)[reference_index])
anchors <- FindIntegrationAnchors(
  object.list = objects,
  normalization.method = "SCT",
  anchor.features = features,
  reduction = "rpca",
  reference = reference_index,
  dims = 1:30,
  k.anchor = 5,
  verbose = TRUE
)
integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT", dims = 1:30, verbose = TRUE)
DefaultAssay(integrated) <- "integrated"
integrated <- RunPCA(integrated, npcs = 50, approx = TRUE, verbose = FALSE)
integrated <- FindNeighbors(integrated, dims = 1:30, verbose = FALSE)
integrated <- FindClusters(
  integrated,
  resolution = c(0.3, 0.6, 1.0),
  algorithm = 1,
  random.seed = 20260726,
  verbose = FALSE
)
primary_column <- grep("_res\\.0\\.6$", colnames(integrated[[]]), value = TRUE)
if (length(primary_column) != 1L) stop("Could not identify resolution 0.6 cluster column")
Idents(integrated) <- integrated[[primary_column, drop = TRUE]]
integrated$seurat_clusters <- as.character(Idents(integrated))
integrated <- RunUMAP(integrated, dims = 1:30, seed.use = 20260726, verbose = FALSE)
integrated <- tryCatch(JoinLayers(integrated, assay = "RNA"), error = function(e) integrated)

saveRDS(integrated, file.path(result_dir, "10samples_integrated_seurat_v5.rds"), compress = FALSE)

embedding <- as.data.frame(Embeddings(integrated, "umap"))
embedding$cell_id <- rownames(embedding)
embedding$sample_id <- integrated$sample_id
embedding$severity <- integrated$severity
embedding$cluster <- integrated$seurat_clusters
write.csv(embedding, gzfile(file.path(result_dir, "cell_embeddings.csv.gz")), row.names = FALSE)

cluster_counts <- as.data.frame(table(
  cluster = integrated$seurat_clusters,
  sample_id = integrated$sample_id,
  severity = integrated$severity
), stringsAsFactors = FALSE)
cluster_counts <- cluster_counts[cluster_counts$Freq > 0, , drop = FALSE]
write.csv(cluster_counts, file.path(result_dir, "cluster_counts_by_sample.csv"), row.names = FALSE)

count_matrix <- xtabs(Freq ~ cluster + sample_id, data = cluster_counts)
entropy <- apply(count_matrix, 1, function(values) {
  probabilities <- values / sum(values)
  probabilities <- probabilities[probabilities > 0]
  -sum(probabilities * log(probabilities)) / log(length(samples))
})
balance <- data.frame(
  cluster = rownames(count_matrix),
  cells = rowSums(count_matrix),
  samples_present = rowSums(count_matrix > 0),
  normalized_sample_entropy = entropy,
  stringsAsFactors = FALSE
)
write.csv(balance, file.path(result_dir, "integration_balance.csv"), row.names = FALSE)

p_cluster <- DimPlot(integrated, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE, raster = TRUE) + NoLegend()
p_sample <- DimPlot(integrated, reduction = "umap", group.by = "sample_id", raster = TRUE)
p_severity <- DimPlot(integrated, reduction = "umap", group.by = "severity", raster = TRUE, cols = c(mild = "#2878B5", severe = "#D64B40"))
ggsave(file.path(result_dir, "umap_clusters.pdf"), p_cluster, width = 9, height = 7)
ggsave(file.path(result_dir, "umap_clusters.png"), p_cluster, width = 9, height = 7, dpi = 170)
ggsave(file.path(result_dir, "umap_samples.pdf"), p_sample, width = 10, height = 7)
ggsave(file.path(result_dir, "umap_severity.pdf"), p_severity, width = 9, height = 7)

methods <- c(
  "# Stage 02: SCTransform and RPCA integration",
  "",
  paste0("- Cells: ", format(ncol(integrated), big.mark = ",")),
  paste0("- Integration features: ", length(features)),
  paste0("- RPCA reference sample: ", names(objects)[reference_index]),
  paste0("- Parallelism: ", sample_workers, " independent SCTransform sample workers; RPCA integration remained one joint task."),
  "- SCTransform was fitted independently within each patient sample using glmGamPoi.",
  "- RNA counts were read directly from each filtered RSEC MEX archive after QC approval.",
  "- RPCA integration was used only for PCA, neighbors, UMAP, and clustering.",
  "- Severity was retained as biological metadata and was not regressed or integrated away.",
  "- Unmodified RNA counts remain in the RNA assay for pseudobulk aggregation.",
  "- Primary clustering resolution: 0.6; sensitivity resolutions: 0.3 and 1.0."
)
writeLines(methods, file.path(result_dir, "integration_report.md"))
writeLines(capture.output(sessionInfo()), file.path(result_dir, "r_sessionInfo.txt"))
