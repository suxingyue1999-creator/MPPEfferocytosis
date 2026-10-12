#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
})

set.seed(20260721)
options(future.globals.maxSize = 40 * 1024^3)
future::plan("sequential")

workspace <- Sys.getenv("MPP_WORKSPACE", "inputs/workspace")
work_root <- file.path(workspace, "工作文件", "23例单细胞重新分析")
input_path <- file.path(
  work_root,
  "results",
  "02_rna_reanalysis",
  "23samples_RNA_reclustered_annotated.rds"
)
out_dir <- file.path(work_root, "results", "02_rna_reanalysis", "subclusters")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

required_meta <- c(
  "sample_id", "severity", "allergy_group", "plastic_group",
  "batch_wave", "celltype_broad_rechecked"
)

sanitize_name <- function(x) {
  gsub("[^A-Za-z0-9_]+", "_", x)
}

join_rna_layers <- function(object) {
  DefaultAssay(object) <- "RNA"
  layer_names <- tryCatch(Layers(object[["RNA"]]), error = function(e) character())
  if (length(layer_names) > 1L) {
    object <- JoinLayers(object, assay = "RNA")
  }
  object
}

build_clean_lineage_object <- function(object, cells) {
  lineage <- subset(object, cells = cells)
  lineage <- join_rna_layers(lineage)
  counts <- LayerData(lineage[["RNA"]], layer = "counts")
  meta <- lineage[[]]
  meta <- meta[colnames(counts), , drop = FALSE]
  clean <- CreateSeuratObject(counts = counts, meta.data = meta)
  DefaultAssay(clean) <- "RNA"
  clean
}

choose_npcs <- function(object, variable_features) {
  max_pcs <- min(
    50L,
    length(variable_features),
    max(5L, ncol(object) - 1L)
  )
  max(5L, max_pcs)
}

choose_dims <- function(object, reduction = "pca") {
  stdev <- Stdev(object, reduction = reduction)
  dims <- sum(is.finite(stdev) & stdev > 0)
  data_layers <- Layers(object[["RNA"]], search = "^data\\.")
  if (length(data_layers)) {
    layer_cells <- vapply(
      data_layers,
      function(layer) ncol(LayerData(object[["RNA"]], layer = layer)),
      integer(1)
    )
    dims <- min(dims, min(layer_cells) - 1L)
  }
  if (dims < 1L) stop("No usable PCs were produced for reduction: ", reduction)
  seq_len(min(30L, dims))
}

pick_resolution <- function(cell_count) {
  if (cell_count < 1500L) return(0.3)
  if (cell_count < 5000L) return(0.5)
  0.7
}

downsample_cells <- function(groups, max_per_group = 500L) {
  pieces <- split(names(groups), groups)
  unlist(lapply(pieces, function(cells) {
    if (length(cells) <= max_per_group) {
      cells
    } else {
      sample(cells, size = max_per_group)
    }
  }), use.names = FALSE)
}

table_to_long <- function(tab, value_name) {
  out <- as.data.frame(tab, stringsAsFactors = FALSE)
  colnames(out)[ncol(out)] <- value_name
  out
}

compute_group_composition <- function(meta, label_col, group_vars) {
  counts_list <- vector("list", length(group_vars))
  fractions_list <- vector("list", length(group_vars))
  for (i in seq_along(group_vars)) {
    group_var <- group_vars[[i]]
    tab <- table(meta[[group_var]], meta[[label_col]])
    count_df <- table_to_long(tab, "cell_count")
    colnames(count_df)[1:2] <- c("group_level", label_col)
    count_df$group_variable <- group_var
    counts_list[[i]] <- count_df[, c("group_variable", "group_level", label_col, "cell_count")]

    frac_df <- count_df
    totals <- ave(frac_df$cell_count, frac_df$group_level, FUN = sum)
    frac_df$cell_fraction <- ifelse(totals > 0, frac_df$cell_count / totals, NA_real_)
    fractions_list[[i]] <- frac_df[, c("group_variable", "group_level", label_col, "cell_count", "cell_fraction")]
  }
  list(
    counts = do.call(rbind, counts_list),
    fractions = do.call(rbind, fractions_list)
  )
}

assign_cluster_labels <- function(meta, cluster_col, module_cols, legacy_col = NULL) {
  cluster_ids <- sort(unique(as.character(meta[[cluster_col]])))
  rows <- lapply(cluster_ids, function(cluster_id) {
    index <- which(as.character(meta[[cluster_col]]) == cluster_id)
    module_labels <- meta$module_majority_call[index]
    module_table <- sort(table(module_labels), decreasing = TRUE)
    module_majority <- names(module_table)[1]
    module_purity <- as.numeric(module_table[1]) / sum(module_table)
    score_margin <- median(meta$module_top_margin[index], na.rm = TRUE)

    legacy_majority <- NA_character_
    legacy_purity <- NA_real_
    if (!is.null(legacy_col) && legacy_col %in% colnames(meta)) {
      legacy_values <- meta[[legacy_col]][index]
      legacy_values <- legacy_values[!is.na(legacy_values) & nzchar(legacy_values)]
      if (length(legacy_values)) {
        legacy_table <- sort(table(legacy_values), decreasing = TRUE)
        legacy_majority <- names(legacy_table)[1]
        legacy_purity <- as.numeric(legacy_table[1]) / sum(legacy_table)
      }
    }

    data.frame(
      cluster_id = cluster_id,
      cells = length(index),
      module_majority = module_majority,
      module_purity = module_purity,
      median_module_margin = score_margin,
      legacy_majority = legacy_majority,
      legacy_purity = legacy_purity,
      final_subtype = module_majority,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

top_markers_per_cluster <- function(markers, n = 10L) {
  if (!nrow(markers)) return(markers)
  pieces <- split(markers, markers$cluster)
  out <- lapply(pieces, function(df) {
    df <- df[order(-df$avg_log2FC, df$p_val_adj, df$gene), , drop = FALSE]
    head(df, n)
  })
  do.call(rbind, out)
}

lineage_configs <- list(
  list(
    lineage_name = "Macrophage_Monocyte",
    broad_labels = c("Macrophage", "Monocyte"),
    legacy_col = "previous_celltype_fine",
    modules = list(
      FABP4_alveolar = c("FABP4", "PPARG", "MARCO", "C1QA", "C1QB", "APOC1", "INHBA"),
      FCN1_monocyte = c("FCN1", "VCAN", "CTSS", "S100A8", "S100A9", "LILRB1", "SAT1"),
      inflammatory = c("IL1B", "CXCL8", "CCL2", "CXCL2", "TNF", "NFKBIA", "EREG"),
      IFN_response = c("ISG15", "IFI6", "IFIT1", "IFIT3", "IFITM3", "MX1", "IFI44L"),
      proliferating = c("MKI67", "TOP2A", "TYMS", "STMN1", "HMGB2", "PCNA")
    )
  ),
  list(
    lineage_name = "T_cells",
    broad_labels = c("T_cells"),
    legacy_col = "previous_celltype_fine",
    modules = list(
      IL7R_naive_like = c("IL7R", "LTB", "TCF7", "LEF1", "MAL", "LTB"),
      cytotoxic = c("NKG7", "CCL5", "CTSW", "PRF1", "GZMB", "GNLY", "KLRD1"),
      Treg_like = c("IL2RA", "FOXP3", "TIGIT", "CTLA4", "IKZF2", "TNFRSF4"),
      IFN_response = c("ISG15", "IFI6", "IFIT1", "IFIT3", "MX1", "OASL"),
      proliferating = c("MKI67", "TOP2A", "TYMS", "STMN1", "HMGB2", "PCNA")
    )
  ),
  list(
    lineage_name = "B_cell",
    broad_labels = c("B_cell"),
    legacy_col = "previous_celltype_fine",
    modules = list(
      naive_B = c("MS4A1", "CD79A", "CD79B", "CD74", "HVCN1", "BANK1", "BLNK"),
      plasma_B = c("MZB1", "JCHAIN", "IGKC", "XBP1", "PRDM1", "SSR4"),
      activated_B = c("CD83", "CD69", "HLA-DRA", "HLA-DPA1", "CD40", "NFKBIA"),
      IFN_response = c("ISG15", "IFI6", "IFIT1", "IFIT3", "MX1", "IFI44L"),
      proliferating = c("MKI67", "TOP2A", "TYMS", "STMN1", "HMGB2", "PCNA")
    )
  ),
  list(
    lineage_name = "Neutrophils",
    broad_labels = c("Neutrophils"),
    legacy_col = "previous_celltype_fine",
    modules = list(
      mature_neutrophil = c("FCGR3B", "CXCR2", "CSF3R", "FPR1", "MNDA", "SIGLEC5"),
      inflammatory = c("S100A8", "S100A9", "CXCL8", "IL1B", "NAMPT", "NFKBIA"),
      immature_neutrophil = c("RETN", "MMP9", "CAMP", "ELANE", "BPI", "DEFA3"),
      IFN_response = c("ISG15", "IFI6", "IFIT1", "IFIT3", "MX1", "OASL"),
      proliferating = c("MKI67", "TOP2A", "TYMS", "STMN1", "HMGB2", "PCNA")
    )
  )
)

lineage_only <- Sys.getenv("LINEAGE_ONLY", unset = "")
if (nzchar(lineage_only)) {
  lineage_configs <- Filter(
    function(cfg) identical(cfg$lineage_name, lineage_only),
    lineage_configs
  )
  if (!length(lineage_configs)) stop("Unknown LINEAGE_ONLY value: ", lineage_only)
  message("Running only requested lineage: ", lineage_only)
}

if (!file.exists(input_path)) {
  stop(
    "Annotated RNA object is missing: ", input_path, "\n",
    "Run scripts/02_cluster_annotate_rna.R first."
  )
}

message("Loading annotated RNA object")
object <- readRDS(input_path)
if (!inherits(object, "Seurat")) stop("Input is not a Seurat object: ", input_path)
DefaultAssay(object) <- "RNA"

missing_meta <- setdiff(required_meta, colnames(object[[]]))
if (length(missing_meta)) {
  stop("Required metadata columns are missing: ", paste(missing_meta, collapse = ", "))
}

# Lineages are rebuilt from raw counts below; retaining normalized/scaled data,
# graphs, and reductions from the 209k-cell parent only increases peak memory.
message("Trimming parent object to RNA counts and metadata")
object <- DietSeurat(
  object,
  assays = "RNA",
  layers = "counts",
  dimreducs = NULL,
  graphs = NULL
)
gc(verbose = FALSE)

overview_rows <- vector("list", length(lineage_configs))

for (cfg_index in seq_along(lineage_configs)) {
  cfg <- lineage_configs[[cfg_index]]
  lineage_name <- cfg$lineage_name
  lineage_stem <- sanitize_name(lineage_name)
  lineage_dir <- file.path(out_dir, lineage_stem)
  dir.create(lineage_dir, recursive = TRUE, showWarnings = FALSE)

  message("Subclustering lineage: ", lineage_name)
  lineage_cells <- rownames(object[[]])[object$celltype_broad_rechecked %in% cfg$broad_labels]
  if (length(lineage_cells) < 200L) {
    warning(lineage_name, ": fewer than 200 cells; skipping")
    overview_rows[[cfg_index]] <- data.frame(
      lineage = lineage_name,
      broad_labels = paste(cfg$broad_labels, collapse = ";"),
      cells = length(lineage_cells),
      samples = length(unique(object$sample_id[match(lineage_cells, colnames(object))])),
      clusters = NA_integer_,
      subtypes = NA_integer_,
      status = "skipped_low_cell_count",
      stringsAsFactors = FALSE
    )
    next
  }

  lineage <- build_clean_lineage_object(object, lineage_cells)
  lineage$lineage_source_broad <- paste(sort(unique(cfg$broad_labels)), collapse = "+")
  lineage$sample_id <- as.character(lineage$sample_id)
  lineage$severity <- as.character(lineage$severity)
  lineage$allergy_group <- as.character(lineage$allergy_group)
  lineage$plastic_group <- as.character(lineage$plastic_group)

  sample_count <- length(unique(lineage$sample_id))
  if (sample_count >= 2L) {
    lineage[["RNA"]] <- split(lineage[["RNA"]], f = lineage$sample_id)
  }

  lineage <- NormalizeData(
    lineage,
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = TRUE
  )
  lineage <- FindVariableFeatures(
    lineage,
    selection.method = "vst",
    nfeatures = 3000,
    verbose = TRUE
  )
  variable_features <- VariableFeatures(lineage)
  if (length(variable_features) < 200L) {
    warning(lineage_name, ": too few variable features after filtering; skipping")
    overview_rows[[cfg_index]] <- data.frame(
      lineage = lineage_name,
      broad_labels = paste(cfg$broad_labels, collapse = ";"),
      cells = ncol(lineage),
      samples = sample_count,
      clusters = NA_integer_,
      subtypes = NA_integer_,
      status = "skipped_low_variable_features",
      stringsAsFactors = FALSE
    )
    next
  }

  regress_vars <- intersect(c("percent.mt"), colnames(lineage[[]]))
  lineage <- ScaleData(
    lineage,
    features = variable_features,
    vars.to.regress = regress_vars,
    verbose = TRUE
  )
  npcs <- choose_npcs(lineage, variable_features)
  lineage <- RunPCA(
    lineage,
    features = variable_features,
    npcs = npcs,
    verbose = TRUE
  )
  dims_use <- choose_dims(lineage, reduction = "pca")

  reduction_name <- "pca"
  if (sample_count >= 2L) {
    data_layers <- Layers(lineage[["RNA"]], search = "^data\\.")
    layer_cells <- vapply(
      data_layers,
      function(layer) ncol(LayerData(lineage[["RNA"]], layer = layer)),
      integer(1)
    )
    min_layer_cells <- min(layer_cells)
    integration_k <- max(2L, min(100L, min_layer_cells - 1L))
    integration_weight_k <- max(2L, min(50L, min_layer_cells - 1L))
    integration_anchor_k <- max(2L, min(5L, min_layer_cells - 1L))
    integration_score_k <- max(2L, min(30L, min_layer_cells - 1L))
    preferred_references <- c("V01", "V18", "V07", "V14")
    reference_indices <- which(data_layers %in% paste0("data.", preferred_references))
    if (!length(reference_indices)) {
      reference_indices <- 1L
    }
    lineage <- IntegrateLayers(
      object = lineage,
      method = RPCAIntegration,
      orig.reduction = "pca",
      new.reduction = "integrated.rpca",
      normalization.method = "LogNormalize",
      dims = dims_use,
      reference = reference_indices,
      k.filter = integration_k,
      k.weight = integration_weight_k,
      k.anchor = integration_anchor_k,
      k.score = integration_score_k,
      verbose = TRUE
    )
    reduction_name <- "integrated.rpca"
  }

  lineage <- FindNeighbors(
    lineage,
    reduction = reduction_name,
    dims = dims_use,
    k.param = 30,
    verbose = TRUE
  )
  resolution <- pick_resolution(ncol(lineage))
  lineage <- FindClusters(
    lineage,
    resolution = resolution,
    algorithm = 4,
    leiden_method = "igraph",
    random.seed = 20260721,
    verbose = TRUE
  )
  cluster_col <- paste0("RNA_snn_res.", resolution)
  lineage$subcluster_id <- lineage[[cluster_col]][, 1]
  Idents(lineage) <- "subcluster_id"
  lineage <- RunUMAP(
    lineage,
    reduction = reduction_name,
    dims = dims_use,
    n.neighbors = 30,
    min.dist = 0.3,
    metric = "cosine",
    seed.use = 20260721,
    reduction.name = "umap.lineage",
    verbose = TRUE
  )

  lineage <- JoinLayers(lineage, assay = "RNA")
  DefaultAssay(lineage) <- "RNA"

  module_panels <- lapply(cfg$modules, function(genes) intersect(unique(genes), rownames(lineage)))
  keep_modules <- lengths(module_panels) >= 2L
  module_panels <- module_panels[keep_modules]
  if (!length(module_panels)) {
    warning(lineage_name, ": all subtype marker modules are unavailable")
    overview_rows[[cfg_index]] <- data.frame(
      lineage = lineage_name,
      broad_labels = paste(cfg$broad_labels, collapse = ";"),
      cells = ncol(lineage),
      samples = sample_count,
      clusters = length(unique(lineage$subcluster_id)),
      subtypes = NA_integer_,
      status = "skipped_no_marker_modules",
      stringsAsFactors = FALSE
    )
    next
  }

  lineage <- AddModuleScore(
    lineage,
    features = module_panels,
    assay = "RNA",
    name = "subtype_module_",
    seed = 20260721,
    search = FALSE
  )
  module_cols <- paste0("subtype_module_", seq_along(module_panels))
  module_score_matrix <- as.matrix(lineage[[]][, module_cols, drop = FALSE])
  colnames(module_score_matrix) <- names(module_panels)
  module_score_order <- t(apply(module_score_matrix, 1, function(x) order(x, decreasing = TRUE)))
  lineage$module_majority_call <- colnames(module_score_matrix)[module_score_order[, 1]]
  lineage$module_top_margin <- apply(module_score_matrix, 1, function(x) {
    ordered <- sort(x, decreasing = TRUE)
    if (length(ordered) >= 2L) ordered[1] - ordered[2] else ordered[1]
  })

  cluster_annotation <- assign_cluster_labels(
    meta = lineage[[]],
    cluster_col = "subcluster_id",
    module_cols = colnames(module_score_matrix),
    legacy_col = cfg$legacy_col
  )
  subtype_map <- setNames(cluster_annotation$final_subtype, cluster_annotation$cluster_id)
  lineage$subtype_label <- unname(subtype_map[as.character(lineage$subcluster_id)])

  module_scores_by_cluster <- aggregate(
    module_score_matrix,
    by = list(subcluster_id = lineage$subcluster_id, subtype_label = lineage$subtype_label),
    FUN = mean
  )
  module_scores_by_cluster <- module_scores_by_cluster[
    order(module_scores_by_cluster$subcluster_id),
    ,
    drop = FALSE
  ]

  marker_cells <- downsample_cells(setNames(as.character(lineage$subcluster_id), Cells(lineage)), max_per_group = 500L)
  marker_object <- subset(lineage, cells = marker_cells)
  Idents(marker_object) <- "subcluster_id"
  markers <- FindAllMarkers(
    marker_object,
    assay = "RNA",
    only.pos = TRUE,
    logfc.threshold = 0.25,
    min.pct = 0.1,
    test.use = "wilcox",
    return.thresh = 0.05,
    verbose = TRUE
  )
  if (nrow(markers)) {
    markers$gene <- rownames(markers)
    rownames(markers) <- NULL
    markers <- markers[, c("gene", setdiff(colnames(markers), "gene"))]
  }
  top_markers <- top_markers_per_cluster(markers, n = 10L)

  lineage_meta <- lineage[[]]
  lineage_meta$cell_id <- rownames(lineage_meta)
  composition_subtype <- compute_group_composition(
    meta = lineage_meta,
    label_col = "subtype_label",
    group_vars = c("sample_id", "severity", "allergy_group", "plastic_group")
  )
  composition_cluster <- compute_group_composition(
    meta = lineage_meta,
    label_col = "subcluster_id",
    group_vars = c("sample_id", "severity", "allergy_group", "plastic_group")
  )

  genes_for_dotplot <- unique(unlist(module_panels))
  genes_for_dotplot <- genes_for_dotplot[genes_for_dotplot %in% rownames(lineage)]
  dotplot <- DotPlot(
    lineage,
    features = genes_for_dotplot,
    group.by = "subtype_label",
    assay = "RNA"
  ) +
    RotatedAxis() +
    ggtitle(paste0(lineage_name, " subtype markers"))

  p1 <- DimPlot(
    lineage,
    reduction = "umap.lineage",
    group.by = "subtype_label",
    raster = TRUE,
    label = TRUE,
    repel = TRUE
  ) + ggtitle(paste0(lineage_name, ": subtype labels"))
  p2 <- DimPlot(
    lineage,
    reduction = "umap.lineage",
    group.by = "subcluster_id",
    raster = TRUE,
    label = TRUE,
    repel = TRUE
  ) + ggtitle(paste0(lineage_name, ": RNA subclusters"))
  p3 <- DimPlot(
    lineage,
    reduction = "umap.lineage",
    group.by = "sample_id",
    raster = TRUE
  ) + ggtitle(paste0(lineage_name, ": samples"))
  p4 <- DimPlot(
    lineage,
    reduction = "umap.lineage",
    group.by = "severity",
    raster = TRUE
  ) + ggtitle(paste0(lineage_name, ": severity"))

  ggsave(
    file.path(lineage_dir, paste0(lineage_stem, "_umap_panels.pdf")),
    p1 + p2 + p3 + p4,
    width = 18,
    height = 12,
    useDingbats = FALSE
  )
  ggsave(
    file.path(lineage_dir, paste0(lineage_stem, "_umap_panels.png")),
    p1 + p2 + p3 + p4,
    width = 18,
    height = 12,
    dpi = 180
  )
  ggsave(
    file.path(lineage_dir, paste0(lineage_stem, "_marker_dotplot.pdf")),
    dotplot,
    width = max(10, length(genes_for_dotplot) * 0.35),
    height = 6,
    useDingbats = FALSE
  )

  writeLines(variable_features, file.path(lineage_dir, paste0(lineage_stem, "_variable_features.txt")))
  write.csv(
    cluster_annotation,
    file.path(lineage_dir, paste0(lineage_stem, "_cluster_annotation_summary.csv")),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
  write.csv(
    module_scores_by_cluster,
    file.path(lineage_dir, paste0(lineage_stem, "_module_scores_by_cluster.csv")),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
  write.csv(
    markers,
    file.path(lineage_dir, paste0(lineage_stem, "_cluster_markers_all.csv")),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
  write.csv(
    top_markers,
    file.path(lineage_dir, paste0(lineage_stem, "_top10_markers_by_cluster.csv")),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
  write.csv(
    composition_subtype$counts,
    file.path(lineage_dir, paste0(lineage_stem, "_subtype_counts_by_group.csv")),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
  write.csv(
    composition_subtype$fractions,
    file.path(lineage_dir, paste0(lineage_stem, "_subtype_fractions_by_group.csv")),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
  write.csv(
    composition_cluster$counts,
    file.path(lineage_dir, paste0(lineage_stem, "_cluster_counts_by_group.csv")),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
  write.csv(
    composition_cluster$fractions,
    file.path(lineage_dir, paste0(lineage_stem, "_cluster_fractions_by_group.csv")),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )

  cell_meta_cols <- intersect(
    c(
      "cell_id", "sample_id", "severity", "allergy_group", "plastic_group",
      "batch_wave", "nCount_RNA", "nFeature_RNA", "percent.mt",
      "previous_celltype_main", "previous_celltype_fine", "cluster_rna",
      "subcluster_id", "module_majority_call", "module_top_margin",
      "subtype_label", "legacy_annotation_usable", "ADT_matrix_usable"
    ),
    colnames(lineage_meta)
  )
  cell_con <- gzfile(file.path(lineage_dir, paste0(lineage_stem, "_cell_metadata.tsv.gz")), open = "wt")
  write.table(
    lineage_meta[, cell_meta_cols, drop = FALSE],
    file = cell_con,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE,
    na = "NA"
  )
  close(cell_con)

  saveRDS(
    lineage,
    file.path(lineage_dir, paste0(lineage_stem, "_subclustered.rds")),
    compress = "gzip"
  )

  subtype_counts <- sort(table(lineage$subtype_label), decreasing = TRUE)
  report_lines <- c(
    paste0("# ", lineage_name, " 亚群重分析"),
    "",
    paste0("- 运行时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste0("- 输入宽谱细胞类型：", paste(cfg$broad_labels, collapse = ", ")),
    paste0("- 细胞数：", format(ncol(lineage), big.mark = ",")),
    paste0("- 样本数：", sample_count),
    paste0("- 可变基因数：", length(variable_features)),
    paste0("- 聚类分辨率：", resolution),
    paste0("- RNA 子群数：", length(unique(lineage$subcluster_id))),
    paste0("- 最终亚群标签数：", length(unique(lineage$subtype_label))),
    "- 方法：从主对象提取该谱系细胞后，仅使用 RNA 原始 counts 重建对象；按 sample_id 分层 LogNormalize，RPCA 仅用于批次整合后的邻居图/聚类/UMAP。",
    "- 注释：以预设标志基因模块得分多数标签定义亚群，旧 `previous_celltype_fine` 只作为参考，不直接继承。",
    "- 解释边界：本步骤仅使用 RNA；当前矩阵中的 ADT/WNN 已另行重建，旧供应商 RDS/WNN 标签不在本步骤复用。",
    "",
    "## 亚群规模",
    ""
  )
  report_lines <- c(
    report_lines,
    paste0("- ", names(subtype_counts), ": ", as.integer(subtype_counts))
  )
  report_lines <- c(
    report_lines,
    "",
    "`*_cluster_annotation_summary.csv`、`*_cluster_markers_all.csv` 和 `*_subtype_fractions_by_group.csv` 分别保留聚类到亚群映射、cluster marker 以及 sample/临床组构成。"
  )
  writeLines(
    report_lines,
    file.path(lineage_dir, paste0(lineage_stem, "_analysis_report.md")),
    useBytes = TRUE
  )

  overview_rows[[cfg_index]] <- data.frame(
    lineage = lineage_name,
    broad_labels = paste(cfg$broad_labels, collapse = ";"),
    cells = ncol(lineage),
    samples = sample_count,
    clusters = length(unique(lineage$subcluster_id)),
    subtypes = length(unique(lineage$subtype_label)),
    status = "completed",
    stringsAsFactors = FALSE
  )
}

overview <- do.call(rbind, overview_rows)
overview_path <- if (nzchar(lineage_only)) {
  file.path(out_dir, paste0("subcluster_overview_", sanitize_name(lineage_only), ".csv"))
} else {
  file.path(out_dir, "subcluster_overview.csv")
}
write.csv(
  overview,
  overview_path,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

writeLines(
  c(
    "# 23例四大谱系亚群重分析",
    "",
    paste0("- 运行时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "- 处理谱系：Macrophage+Monocyte、T_cells、B_cell、Neutrophils。",
    "- 主输入：`results/02_rna_reanalysis/23samples_RNA_reclustered_annotated.rds`。",
    "- 核心约束：本步骤仅使用 RNA；V20 保留于 RNA 分析，当前矩阵中的 ADT/WNN 由独立脚本重建。",
    "",
    paste0("`", basename(overview_path), "` 汇总本次运行的谱系细胞数、样本数、子群数和执行状态。")
  ),
  file.path(out_dir, "subcluster_overview.md"),
  useBytes = TRUE
)

session_con <- file(file.path(out_dir, "sessionInfo.txt"), open = "wt")
writeLines(capture.output(sessionInfo()), session_con)
close(session_con)

message("Lineage subclustering script completed")
