#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(hdf5r)
})

samples <- c("D01", "D02", "D03", "D04", "D05", "D06", "D07", "D08", "D09", "D10")
input_root <- Sys.getenv("BD_INPUT_ROOT", "/inputs")
analysis_root <- Sys.getenv("BD_ANALYSIS_ROOT", "/analysis")
audit_dir <- file.path(analysis_root, "results", "00_audit")
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

mex_audit <- read.csv(file.path(audit_dir, "mex_h5ad_audit.csv"), stringsAsFactors = FALSE)
as_bool <- function(value) {
  if (is.logical(value)) return(value)
  tolower(trimws(as.character(value))) == "true"
}
rows <- vector("list", length(samples))
h5_rows <- vector("list", length(samples))

for (index in seq_along(samples)) {
  sample_id <- samples[[index]]
  path <- file.path(input_root, sample_id, paste0(sample_id, "_Seurat.rds"))
  message("Auditing Seurat RDS: ", sample_id)
  object <- readRDS(path)
  if (!inherits(object, "Seurat")) stop(sample_id, ": RDS is not a Seurat object")
  default_assay <- DefaultAssay(object)
  counts <- GetAssayData(object, assay = default_assay, slot = "counts")
  expected <- mex_audit[mex_audit$sample_id == sample_id, , drop = FALSE]
  if (nrow(expected) != 1L) stop(sample_id, ": MEX audit row missing")
  rows[[index]] <- data.frame(
    sample_id = sample_id,
    object_class = paste(class(object), collapse = "|"),
    object_version = as.character(object@version),
    default_assay = default_assay,
    assays = paste(Assays(object), collapse = "|"),
    cells = ncol(object),
    features = nrow(object),
    count_rows = nrow(counts),
    count_columns = ncol(counts),
    cells_unique = !anyDuplicated(Cells(object)),
    features_unique = !anyDuplicated(rownames(object)),
    cells_equal_filtered_mex = ncol(object) == expected$filtered_cells,
    features_equal_filtered_mex = nrow(object) == expected$features,
    stringsAsFactors = FALSE
  )
  rm(object, counts)
  gc(verbose = FALSE)

  h5_path <- file.path(input_root, sample_id, paste0(sample_id, ".h5ad"))
  message("Auditing H5AD structure: ", sample_id)
  handle <- H5File$new(h5_path, mode = "r")
  x <- handle[["X"]]
  shape <- as.integer(h5attr(x, "shape"))
  obs <- handle[["obs"]]
  var <- handle[["var"]]
  obs_index_name <- as.character(h5attr(obs, "_index"))
  var_index_name <- as.character(h5attr(var, "_index"))
  obs_names <- as.character(obs[[obs_index_name]][])
  var_names <- as.character(var[[var_index_name]][])
  h5_rows[[index]] <- data.frame(
    sample_id = sample_id,
    n_obs = shape[[1]],
    n_vars = shape[[2]],
    x_rows = shape[[1]],
    x_columns = shape[[2]],
    obs_names_unique = !anyDuplicated(obs_names),
    var_names_unique = !anyDuplicated(var_names),
    obs_columns = paste(names(obs), collapse = "|"),
    var_columns = paste(names(var), collapse = "|"),
    obsm_keys = paste(names(handle[["obsm"]]), collapse = "|"),
    uns_keys = paste(names(handle[["uns"]]), collapse = "|"),
    cells_equal_filtered_mex = shape[[1]] == expected$filtered_cells,
    features_equal_filtered_mex = shape[[2]] == expected$features,
    stringsAsFactors = FALSE
  )
  handle$close_all()
  rm(handle, x, obs, var, obs_names, var_names)
  gc(verbose = FALSE)
}

audit <- do.call(rbind, rows)
h5_audit <- do.call(rbind, h5_rows)
write.csv(audit, file.path(audit_dir, "seurat_rds_audit.csv"), row.names = FALSE)
write.csv(h5_audit, file.path(audit_dir, "h5ad_schema_audit.csv"), row.names = FALSE)
mex_audit$h5ad_cells <- h5_audit$n_obs[match(mex_audit$sample_id, h5_audit$sample_id)]
mex_audit$h5ad_features <- h5_audit$n_vars[match(mex_audit$sample_id, h5_audit$sample_id)]
mex_audit$all_cross_format_checks_pass <- as_bool(mex_audit$all_cross_format_checks_pass) &
  h5_audit$cells_equal_filtered_mex[match(mex_audit$sample_id, h5_audit$sample_id)] &
  h5_audit$features_equal_filtered_mex[match(mex_audit$sample_id, h5_audit$sample_id)] &
  h5_audit$obs_names_unique[match(mex_audit$sample_id, h5_audit$sample_id)] &
  h5_audit$var_names_unique[match(mex_audit$sample_id, h5_audit$sample_id)]
write.csv(mex_audit, file.path(audit_dir, "mex_h5ad_audit.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(audit_dir, "r_sessionInfo.txt"))

checks <- c(
  audit$cells_unique,
  audit$features_unique,
  audit$cells_equal_filtered_mex,
  audit$features_equal_filtered_mex
  , h5_audit$obs_names_unique
  , h5_audit$var_names_unique
  , h5_audit$cells_equal_filtered_mex
  , h5_audit$features_equal_filtered_mex
)
if (!all(checks)) stop("Seurat RDS/H5AD cross-format audit failed")
message("Seurat RDS and H5AD audits completed")
