#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(data.table)
  library(future)
})

set.seed(20260820)
options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 20 * 1024^3)

workspace <- Sys.getenv("MPP_WORKSPACE", "inputs/workspace")
work_root <- file.path(workspace, "工作文件")
analysis_root <- file.path(work_root, "23例单细胞重新分析")
input_path <- file.path(analysis_root, "results", "02_rna_reanalysis", "23samples_RNA_reclustered_annotated.rds")
out_dir <- file.path(analysis_root, "results", "06_cellchat_by_sample")
object_dir <- file.path(out_dir, "objects")
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_path)) stop("Missing current RNA object: ", input_path)

message("Loading current RNA object")
object <- readRDS(input_path)
meta <- as.data.table(object@meta.data, keep.rownames = "cell_id")
required <- c("sample_id", "severity", "plastic_group", "allergy_group", "batch_wave", "celltype_broad_rechecked")
missing_cols <- setdiff(required, names(meta))
if (length(missing_cols)) stop("Missing metadata columns: ", paste(missing_cols, collapse = ", "))

meta[, celltype := fifelse(
  celltype_broad_rechecked %chin% c("B_cell", "Macrophage", "Monocyte", "Neutrophils", "T_cells", "Epithelial_cells"),
  celltype_broad_rechecked,
  NA_character_
)]
meta <- meta[!is.na(sample_id) & !is.na(celltype) & !is.na(severity)]
meta[, severity_group := fifelse(severity %chin% c("重症", "severe"), "Severe", "Mild")]
meta[, plastic_group2 := fifelse(plastic_group %chin% c("是", "塑形", "S", "GroupS"), "Plastic", "NonPlastic")]
meta[, allergy_group2 := fifelse(allergy_group %chin% c("是", "过敏", "Allergy"), "Allergy", "NonAllergy")]

sample_info <- unique(meta[, .(sample_id, severity_group, plastic_group2, allergy_group2, batch_wave)])
setorder(sample_info, sample_id)
if (anyDuplicated(sample_info$sample_id)) stop("sample_id has duplicated clinical rows")
fwrite(sample_info, file.path(out_dir, "sample_info.csv"))

counts <- meta[, .(cells = .N), by = .(sample_id, celltype)]
fwrite(counts, file.path(out_dir, "sample_celltype_counts.csv"))

db <- CellChatDB.human
db_genes <- unique(c(db$interaction$ligand, db$interaction$receptor))
db_genes <- db_genes[!grepl("_", db_genes, fixed = TRUE)]
db_genes <- intersect(db_genes, rownames(object))
if (length(db_genes) < 100) stop("Too few CellChat database genes in current RNA object: ", length(db_genes))

data_norm <- GetAssayData(object, assay = "RNA", layer = "data")
data_norm <- data_norm[db_genes, , drop = FALSE]
rm(object)
gc()

run_one <- function(sample_name) {
  path <- file.path(object_dir, paste0(sample_name, ".rds"))
  if (file.exists(path)) {
    message("Reusing ", sample_name)
    return(readRDS(path))
  }
  current <- meta[sample_id == sample_name]
  cells <- current$cell_id
  expr <- data_norm[, cells, drop = FALSE]
  meta_input <- data.frame(
    celltype = factor(current$celltype),
    row.names = cells,
    check.names = FALSE
  )
  cc <- createCellChat(object = expr, meta = meta_input, group.by = "celltype")
  cc@DB <- db
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(
    cc,
    type = "truncatedMean",
    trim = 0.1,
    raw.use = TRUE,
    population.size = FALSE
  )
  cc <- filterCommunication(cc, min.cells = 20)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  saveRDS(cc, path, compress = FALSE)
  cc
}

sample_names <- sample_info$sample_id
log_path <- file.path(out_dir, "run.log")
log_con <- file(log_path, open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")
on.exit({
  sink(type = "message")
  sink(type = "output")
  close(log_con)
}, add = TRUE)

message("CellChat version: ", as.character(packageVersion("CellChat")))
message("Samples: ", length(sample_names), "; database genes retained: ", length(db_genes))
future::plan(sequential)

communications <- vector("list", length(sample_names))
for (i in seq_along(sample_names)) {
  s <- sample_names[[i]]
  message(sprintf("[%d/%d] %s", i, length(sample_names), s))
  cc <- run_one(s)
  x <- tryCatch(as.data.table(subsetCommunication(cc)), error = function(e) data.table())
  if (nrow(x)) {
    x[, sample_id := s]
    communications[[i]] <- x
  }
  rm(cc)
  gc()
}

communication <- rbindlist(communications, use.names = TRUE, fill = TRUE)
if (!nrow(communication)) stop("No communication records returned")
communication <- merge(communication, sample_info, by = "sample_id", all.x = TRUE)
fwrite(communication, file.path(out_dir, "communication_by_sample.csv"))

summarise_feature <- function(dt, keys) {
  x <- dt[, .(
    probability_sum = sum(prob, na.rm = TRUE),
    interaction_count = .N
  ), by = c("sample_id", keys)]
  x <- merge(x, sample_info, by = "sample_id", all.x = TRUE)
  x
}

run_group_test <- function(x, grouping, feature_keys) {
  if (!nrow(x)) return(data.table())
  x[, group := get(grouping)]
  out <- x[, {
    y <- probability_sum
    g <- group
    p <- tryCatch(wilcox.test(y ~ g, exact = FALSE)$p.value, error = function(e) NA_real_)
    data.table(
      n_group1 = sum(g == unique(g)[1]),
      n_group2 = sum(g == unique(g)[2]),
      mean_group1 = mean(y[g == unique(g)[1]], na.rm = TRUE),
      mean_group2 = mean(y[g == unique(g)[2]], na.rm = TRUE),
      log2_ratio = log2((mean(y[g == unique(g)[2]], na.rm = TRUE) + 1e-12) /
                        (mean(y[g == unique(g)[1]], na.rm = TRUE) + 1e-12)),
      p_value = p
    )
  }, by = feature_keys]
  out[, fdr := p.adjust(p_value, method = "BH")]
  out[]
}

pathway <- summarise_feature(communication, "pathway_name")
fwrite(pathway, file.path(out_dir, "pathway_by_sample.csv"))
source_target <- summarise_feature(communication, c("source", "target"))
fwrite(source_target, file.path(out_dir, "source_target_by_sample.csv"))

pathway_severity <- run_group_test(pathway, "severity_group", c("pathway_name"))
pathway_plastic <- run_group_test(pathway, "plastic_group2", c("pathway_name"))
pathway_allergy <- run_group_test(pathway, "allergy_group2", c("pathway_name"))
fwrite(pathway_severity, file.path(out_dir, "pathway_severity_comparison.csv"))
fwrite(pathway_plastic, file.path(out_dir, "pathway_plasticity_comparison.csv"))
fwrite(pathway_allergy, file.path(out_dir, "pathway_allergy_comparison.csv"))

st_severity <- run_group_test(source_target, "severity_group", c("source", "target"))
st_plastic <- run_group_test(source_target, "plastic_group2", c("source", "target"))
st_allergy <- run_group_test(source_target, "allergy_group2", c("source", "target"))
fwrite(st_severity, file.path(out_dir, "source_target_severity_comparison.csv"))
fwrite(st_plastic, file.path(out_dir, "source_target_plasticity_comparison.csv"))
fwrite(st_allergy, file.path(out_dir, "source_target_allergy_comparison.csv"))

audit <- data.table(
  item = c("samples", "cells", "cellchat_database_genes", "communication_records", "sample_objects", "severity_pathways_fdr05", "plasticity_pathways_fdr05", "allergy_pathways_fdr05"),
  value = c(
    length(sample_names), nrow(meta), length(db_genes), nrow(communication),
    sum(file.exists(file.path(object_dir, paste0(sample_names, ".rds")))),
    sum(pathway_severity$fdr < 0.05, na.rm = TRUE),
    sum(pathway_plastic$fdr < 0.05, na.rm = TRUE),
    sum(pathway_allergy$fdr < 0.05, na.rm = TRUE)
  )
)
fwrite(audit, file.path(out_dir, "run_audit.csv"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
message("Completed sample-level CellChat")
