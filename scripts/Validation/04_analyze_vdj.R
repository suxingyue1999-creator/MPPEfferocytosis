#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

options(stringsAsFactors = FALSE)

workspace <- Sys.getenv("MPP_WORKSPACE", "inputs/workspace")
work_root <- file.path(workspace, "工作文件", "23例单细胞重新分析")
out_dir <- file.path(work_root, "results", "04_vdj")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

clinical_path <- file.path(workspace, "工作文件", "单细胞23例样本临床信息.csv")
bcr_base <- file.path(
  workspace,
  "scRNA-seq_data260717",
  "WTD2024-04-006_BCR",
  "BCR",
  "0_Rhapsody"
)
tcr_base <- file.path(
  workspace,
  "scRNA-seq_data260717",
  "WTD2024-04-006_TCR",
  "TCR",
  "0_Rhapsody"
)

metadata_candidates <- c(
  file.path(work_root, "results", "03_subtype", "cell_metadata_subtypes.tsv.gz"),
  file.path(work_root, "results", "02_rna_reanalysis", "cell_metadata_clustered.tsv.gz"),
  file.path(work_root, "results", "02_rna_reanalysis", "cell_metadata_rna_rechecked.tsv.gz"),
  file.path(work_root, "results", "01_qc", "cell_metadata_after_qc.tsv.gz")
)

read_text_table <- function(path, ...) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    read.delim(gzfile(path), sep = "\t", check.names = FALSE, ...)
  } else {
    read.delim(path, sep = "\t", check.names = FALSE, ...)
  }
}

write_tsv_gz <- function(data, path) {
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  write.table(
    data,
    file = con,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE,
    na = ""
  )
}

pick_first_existing <- function(paths) {
  matched <- paths[file.exists(paths)]
  if (length(matched) == 0L) return(NA_character_)
  matched[[1]]
}

file_has_null_bytes <- function(path, bytes = 4096L) {
  if (!file.exists(path)) return(NA)
  con <- file(path, open = "rb")
  on.exit(close(con), add = TRUE)
  raw_block <- readBin(con, what = "raw", n = bytes)
  any(raw_block == as.raw(0))
}

file_header_ok <- function(path) {
  if (!file.exists(path)) return(FALSE)
  if (isTRUE(file_has_null_bytes(path))) return(FALSE)
  ok <- FALSE
  try({
    cols <- names(read.delim(
      path,
      sep = "\t",
      nrows = 0,
      quote = "",
      comment.char = "",
      stringsAsFactors = FALSE,
      fill = TRUE,
      check.names = FALSE
    ))
    ok <- all(c("cell_id", "locus", "productive", "cdr3_aa") %in% cols)
  }, silent = TRUE)
  ok
}

read_airr_subset <- function(path, sample_id) {
  header <- read.delim(
    path,
    sep = "\t",
    nrows = 0,
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    fill = TRUE,
    check.names = FALSE
  )
  needed <- c(
    "cell_id", "locus", "sequence_id", "consensus_count", "duplicate_count",
    "productive", "complete_vdj", "v_call", "d_call", "j_call", "c_call",
    "cdr3", "cdr3_aa"
  )
  available <- intersect(needed, names(header))
  col_classes <- rep("NULL", length(names(header)))
  names(col_classes) <- names(header)
  col_classes[available] <- "character"
  df <- read.delim(
    path,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    fill = TRUE,
    stringsAsFactors = FALSE,
    colClasses = col_classes,
    check.names = FALSE
  )
  df$sample_id <- sample_id
  df$cell_id <- as.character(df$cell_id)
  df$raw_barcode <- paste0("X", df$cell_id)
  df$productive <- toupper(df$productive) %in% c("TRUE", "T", "1")
  df$complete_vdj <- toupper(df$complete_vdj) %in% c("TRUE", "T", "1")
  df$consensus_count <- suppressWarnings(as.numeric(df$consensus_count))
  df$duplicate_count <- suppressWarnings(as.numeric(df$duplicate_count))
  df$v_gene <- sub("\\*.*$", "", ifelse(is.na(df$v_call), "", df$v_call))
  df$d_gene <- sub("\\*.*$", "", ifelse(is.na(df$d_call), "", df$d_call))
  df$j_gene <- sub("\\*.*$", "", ifelse(is.na(df$j_call), "", df$j_call))
  df$c_gene <- sub("\\*.*$", "", ifelse(is.na(df$c_call), "", df$c_call))
  df$receptor <- ifelse(
    grepl("^IG", df$locus), "BCR",
    ifelse(grepl("^TR", df$locus), "TCR", "OTHER")
  )
  df$cdr3_aa <- ifelse(is.na(df$cdr3_aa), "", trimws(df$cdr3_aa))
  df$cdr3_aa_length <- ifelse(nzchar(df$cdr3_aa), nchar(df$cdr3_aa), NA_integer_)
  df
}

top_per_locus <- function(df) {
  if (nrow(df) == 0L) return(df[0, , drop = FALSE])
  ord <- order(
    df$locus,
    -ifelse(is.na(df$consensus_count), -1, df$consensus_count),
    -ifelse(is.na(df$duplicate_count), -1, df$duplicate_count),
    df$sequence_id
  )
  df <- df[ord, , drop = FALSE]
  df[!duplicated(df$locus), , drop = FALSE]
}

collapse_chain_signature <- function(df) {
  if (nrow(df) == 0L) return(NA_character_)
  sig <- paste(df$locus, df$v_gene, df$j_gene, df$cdr3_aa, sep = "|")
  paste(sig[order(df$locus)], collapse = ";")
}

collapse_field_by_locus <- function(df, field) {
  if (nrow(df) == 0L) return(NA_character_)
  vals <- paste(df$locus, df[[field]], sep = ":")
  paste(vals[order(df$locus)], collapse = ";")
}

expansion_category <- function(size_vec) {
  ifelse(
    is.na(size_vec), NA_character_,
    ifelse(
      size_vec <= 1, "singleton",
      ifelse(size_vec <= 4, "small_2_4",
      ifelse(size_vec <= 9, "medium_5_9", "large_ge10"))
    )
  )
}

calc_diversity <- function(counts) {
  counts <- counts[counts > 0]
  total <- sum(counts)
  richness <- length(counts)
  if (total == 0 || richness == 0) {
    return(data.frame(
      clone_cells = 0,
      unique_clones = 0,
      shannon = NA_real_,
      inv_simpson = NA_real_,
      pielou_evenness = NA_real_,
      clonality = NA_real_,
      max_clone_fraction = NA_real_
    ))
  }
  p <- counts / total
  shannon <- -sum(p * log(p))
  inv_simpson <- 1 / sum(p^2)
  evenness <- if (richness > 1) shannon / log(richness) else 1
  data.frame(
    clone_cells = total,
    unique_clones = richness,
    shannon = shannon,
    inv_simpson = inv_simpson,
    pielou_evenness = evenness,
    clonality = 1 - evenness,
    max_clone_fraction = max(p)
  )
}

clinical <- read.csv(clinical_path, check.names = FALSE)
clinical$sample_id <- clinical$sample_code

sample_ids <- clinical$sample_id
file_audit <- do.call(rbind, lapply(sample_ids, function(sample_id) {
  bcr_path <- file.path(bcr_base, sample_id, paste0(sample_id, "_VDJ_Dominant_Contigs_AIRR.tsv"))
  tcr_path <- file.path(tcr_base, sample_id, paste0(sample_id, "_VDJ_Dominant_Contigs_AIRR.tsv"))
  bcr_exists <- file.exists(bcr_path)
  tcr_exists <- file.exists(tcr_path)
  bcr_md5 <- if (bcr_exists) unname(tools::md5sum(bcr_path)) else NA_character_
  tcr_md5 <- if (tcr_exists) unname(tools::md5sum(tcr_path)) else NA_character_
  bcr_readable <- if (bcr_exists) file_header_ok(bcr_path) else FALSE
  tcr_readable <- if (tcr_exists) file_header_ok(tcr_path) else FALSE
  identical_pair <- isTRUE(bcr_exists && tcr_exists && identical(bcr_md5, tcr_md5))
  chosen_source <- if (tcr_readable) "TCR_dir" else if (bcr_readable) "BCR_dir" else NA_character_
  chosen_path <- if (tcr_readable) tcr_path else if (bcr_readable) bcr_path else NA_character_
  issue <- if (!bcr_exists || !tcr_exists) {
    "missing_file"
  } else if (!bcr_readable && tcr_readable) {
    "BCR_unreadable_use_TCR"
  } else if (!tcr_readable && bcr_readable) {
    "TCR_unreadable_use_BCR"
  } else if (identical_pair) {
    "duplicate_pair_identical"
  } else if (bcr_readable && tcr_readable) {
    "pair_content_mismatch"
  } else {
    "both_unreadable"
  }
  data.frame(
    sample_id = sample_id,
    bcr_path = bcr_path,
    tcr_path = tcr_path,
    bcr_exists = bcr_exists,
    tcr_exists = tcr_exists,
    bcr_readable = bcr_readable,
    tcr_readable = tcr_readable,
    bcr_md5 = bcr_md5,
    tcr_md5 = tcr_md5,
    identical_pair = identical_pair,
    chosen_source = chosen_source,
    chosen_path = chosen_path,
    audit_issue = issue,
    stringsAsFactors = FALSE
  )
}))
write.csv(file_audit, file.path(out_dir, "vdj_file_audit.csv"), row.names = FALSE)

if (any(is.na(file_audit$chosen_path))) {
  stop("At least one sample lacks a readable AIRR file")
}

airr_list <- lapply(seq_len(nrow(file_audit)), function(i) {
  read_airr_subset(file_audit$chosen_path[i], file_audit$sample_id[i])
})
airr <- do.call(rbind, airr_list)
airr <- airr[airr$receptor %in% c("BCR", "TCR"), , drop = FALSE]
airr$sample_id <- factor(airr$sample_id, levels = sample_ids)
airr$sample_id <- as.character(airr$sample_id)
write_tsv_gz(airr, file.path(out_dir, "vdj_contigs_selected.tsv.gz"))

contig_summary <- aggregate(
  cbind(
    total_contigs = rep(1L, nrow(airr)),
    productive_contigs = as.integer(airr$productive),
    complete_productive_contigs = as.integer(airr$productive & airr$complete_vdj)
  ),
  by = list(sample_id = airr$sample_id, receptor = airr$receptor, locus = airr$locus),
  FUN = sum
)
write.csv(contig_summary, file.path(out_dir, "contig_summary_by_sample_receptor_locus.csv"), row.names = FALSE)

productive_contigs <- airr[airr$productive & nzchar(airr$cdr3_aa), , drop = FALSE]

cell_receptor_groups <- split(
  productive_contigs,
  paste(productive_contigs$sample_id, productive_contigs$raw_barcode, productive_contigs$receptor, sep = "||"),
  drop = TRUE
)

cell_receptor_summary <- do.call(rbind, lapply(cell_receptor_groups, function(df) {
  top_df <- top_per_locus(df)
  receptor <- unique(df$receptor)[1]
  loci <- sort(unique(df$locus))
  paired <- if (receptor == "BCR") {
    any(loci == "IGH") && any(loci %in% c("IGK", "IGL"))
  } else {
    any(loci == "TRA") && any(loci == "TRB")
  }
  data.frame(
    sample_id = unique(df$sample_id)[1],
    raw_barcode = unique(df$raw_barcode)[1],
    receptor = receptor,
    productive_contigs = nrow(df),
    productive_loci = paste(loci, collapse = ";"),
    productive_locus_count = length(loci),
    paired_chain = paired,
    clone_signature = collapse_chain_signature(top_df),
    cdr3_aa_by_locus = collapse_field_by_locus(top_df, "cdr3_aa"),
    v_gene_by_locus = collapse_field_by_locus(top_df, "v_gene"),
    j_gene_by_locus = collapse_field_by_locus(top_df, "j_gene"),
    c_gene_by_locus = collapse_field_by_locus(top_df, "c_gene"),
    max_consensus_count = max(top_df$consensus_count, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))

if (is.null(cell_receptor_summary) || nrow(cell_receptor_summary) == 0L) {
  stop("No productive receptor contigs were parsed from the AIRR inputs")
}

cell_receptor_summary$clone_id_within_sample <- paste(
  cell_receptor_summary$sample_id,
  cell_receptor_summary$receptor,
  cell_receptor_summary$clone_signature,
  sep = "::"
)
clone_sizes <- as.data.frame(table(cell_receptor_summary$clone_id_within_sample), stringsAsFactors = FALSE)
names(clone_sizes) <- c("clone_id_within_sample", "clone_size")
cell_receptor_summary <- merge(
  cell_receptor_summary,
  clone_sizes,
  by = "clone_id_within_sample",
  all.x = TRUE,
  sort = FALSE
)
cell_receptor_summary$expansion_category <- expansion_category(cell_receptor_summary$clone_size)
cell_receptor_summary$sample_clone_rank <- ave(
  -cell_receptor_summary$clone_size,
  interaction(cell_receptor_summary$sample_id, cell_receptor_summary$receptor, drop = TRUE),
  FUN = function(x) rank(x, ties.method = "min")
)
write_tsv_gz(cell_receptor_summary, file.path(out_dir, "cell_receptor_summary.tsv.gz"))

wide_by_receptor <- split(cell_receptor_summary, cell_receptor_summary$receptor)
rename_with_prefix <- function(df, prefix) {
  keep <- setdiff(names(df), c("sample_id", "raw_barcode", "receptor"))
  names(df)[match(keep, names(df))] <- paste0(prefix, keep)
  df
}

wide_cell <- Reduce(function(x, y) merge(x, y, by = c("sample_id", "raw_barcode"), all = TRUE), lapply(names(wide_by_receptor), function(receptor) {
  prefix <- paste0(tolower(receptor), "_")
  rename_with_prefix(wide_by_receptor[[receptor]][, setdiff(names(wide_by_receptor[[receptor]]), "receptor"), drop = FALSE], prefix)
}))

metadata_path <- pick_first_existing(metadata_candidates)
metadata_source <- if (is.na(metadata_path)) "none" else basename(metadata_path)
if (!is.na(metadata_path)) {
  metadata_all <- read_text_table(metadata_path, stringsAsFactors = FALSE)
  keep_cols <- intersect(
    c(
      "sample_id", "raw_barcode", "cell_id", "severity",
      "allergy_group", "plastic_group", "previous_celltype_main",
      "previous_celltype_fine", "previous_seurat_cluster",
      "cluster_rna", "celltype_broad_rechecked", "RNA_analysis_usable",
      "legacy_annotation_usable"
    ),
    names(metadata_all)
  )
  metadata_use <- metadata_all[, keep_cols, drop = FALSE]
  wide_cell <- merge(wide_cell, metadata_use, by = c("sample_id", "raw_barcode"), all.x = TRUE, sort = FALSE)
}

clinical_keep <- intersect(
  c("sample_id", "severity", "allergy_group", "plastic_group"),
  names(clinical)
)
clinical_use <- unique(clinical[, clinical_keep, drop = FALSE])
wide_cell <- merge(wide_cell, clinical_use, by = "sample_id", all.x = TRUE, sort = FALSE, suffixes = c("", "_clinical"))
for (field in c("severity", "allergy_group", "plastic_group")) {
  clinical_field <- paste0(field, "_clinical")
  if (field %in% names(wide_cell) && clinical_field %in% names(wide_cell)) {
    needs_fill <- is.na(wide_cell[[field]]) | wide_cell[[field]] == ""
    wide_cell[[field]][needs_fill] <- wide_cell[[clinical_field]][needs_fill]
    wide_cell[[clinical_field]] <- NULL
  }
}
wide_cell$rna_metadata_found <- if ("cell_id" %in% names(wide_cell)) !is.na(wide_cell$cell_id) else FALSE
write_tsv_gz(wide_cell, file.path(out_dir, "vdj_cell_level_with_metadata.tsv.gz"))

sample_receptor_summary <- do.call(rbind, lapply(split(cell_receptor_summary, interaction(cell_receptor_summary$sample_id, cell_receptor_summary$receptor, drop = TRUE)), function(df) {
  data.frame(
    sample_id = unique(df$sample_id)[1],
    receptor = unique(df$receptor)[1],
    cells_with_productive_receptor = nrow(df),
    paired_cells = sum(df$paired_chain),
    paired_fraction = mean(df$paired_chain),
    unique_clones = length(unique(df$clone_id_within_sample)),
    expanded_clone_cells_ge2 = sum(df$clone_size >= 2),
    expanded_clone_cells_ge5 = sum(df$clone_size >= 5),
    max_clone_size = max(df$clone_size),
    max_clone_fraction = max(df$clone_size) / nrow(df),
    stringsAsFactors = FALSE
  )
}))

diversity_summary <- do.call(rbind, lapply(split(cell_receptor_summary, interaction(cell_receptor_summary$sample_id, cell_receptor_summary$receptor, drop = TRUE)), function(df) {
  counts <- table(df$clone_id_within_sample)
  out <- calc_diversity(as.numeric(counts))
  out$sample_id <- unique(df$sample_id)[1]
  out$receptor <- unique(df$receptor)[1]
  out
}))
names(diversity_summary)[names(diversity_summary) == "unique_clones"] <- "unique_clones_from_diversity"
names(diversity_summary)[names(diversity_summary) == "max_clone_fraction"] <- "max_clone_fraction_from_diversity"
sample_receptor_summary <- merge(
  sample_receptor_summary,
  diversity_summary,
  by = c("sample_id", "receptor"),
  all.x = TRUE,
  sort = FALSE
)

if ("rna_metadata_found" %in% names(wide_cell)) {
  link_summary <- do.call(rbind, lapply(split(wide_cell, wide_cell$sample_id), function(df) {
    data.frame(
      sample_id = unique(df$sample_id)[1],
      cells_with_any_vdj = nrow(df),
      cells_linked_to_rna_metadata = sum(df$rna_metadata_found, na.rm = TRUE),
      linkage_fraction = mean(df$rna_metadata_found, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  sample_receptor_summary <- merge(sample_receptor_summary, link_summary, by = "sample_id", all.x = TRUE, sort = FALSE)
  write.csv(link_summary, file.path(out_dir, "vdj_rna_linkage_by_sample.csv"), row.names = FALSE)
}
write.csv(sample_receptor_summary, file.path(out_dir, "sample_receptor_summary.csv"), row.names = FALSE)

clone_lookup <- unique(cell_receptor_summary[, c(
  "sample_id", "receptor", "clone_id_within_sample", "clone_signature",
  "clone_size", "expansion_category"
)])
top_clones <- clone_lookup[order(clone_lookup$receptor, -clone_lookup$clone_size, clone_lookup$sample_id), , drop = FALSE]
top_clones <- top_clones[top_clones$clone_size >= 2, , drop = FALSE]
if (nrow(top_clones) > 0L) {
  top_clones <- top_clones[seq_len(min(200L, nrow(top_clones))), , drop = FALSE]
}

celltype_field <- if ("celltype_broad_rechecked" %in% names(wide_cell)) {
  "celltype_broad_rechecked"
} else if ("previous_celltype_main" %in% names(wide_cell)) {
  "previous_celltype_main"
} else {
  NA_character_
}

if (!is.na(celltype_field)) {
  celltype_map <- wide_cell[, c("sample_id", "raw_barcode", celltype_field), drop = FALSE]
  names(celltype_map)[3] <- "celltype_label"
  clone_celltype <- merge(
    cell_receptor_summary,
    celltype_map,
    by = c("sample_id", "raw_barcode"),
    all.x = TRUE,
    sort = FALSE
  )
  clone_celltype$celltype_label[is.na(clone_celltype$celltype_label) | clone_celltype$celltype_label == ""] <- "Unknown"

  clone_celltype_summary <- aggregate(
    cbind(
      cells = rep(1L, nrow(clone_celltype)),
      expanded_clone_cells_ge2 = as.integer(clone_celltype$clone_size >= 2),
      expanded_clone_cells_ge5 = as.integer(clone_celltype$clone_size >= 5)
    ),
    by = list(
      sample_id = clone_celltype$sample_id,
      receptor = clone_celltype$receptor,
      celltype_label = clone_celltype$celltype_label
    ),
    FUN = sum
  )
  write.csv(clone_celltype_summary, file.path(out_dir, "clone_relation_by_celltype.csv"), row.names = FALSE)

  if (nrow(top_clones) > 0L) {
    top_clone_celltypes <- merge(
      top_clones,
      clone_celltype[, c("sample_id", "receptor", "clone_id_within_sample", "celltype_label"), drop = FALSE],
      by = c("sample_id", "receptor", "clone_id_within_sample"),
      all.x = TRUE,
      sort = FALSE
    )
    top_clone_celltypes <- aggregate(
      rep(1L, nrow(top_clone_celltypes)),
      by = list(
        sample_id = top_clone_celltypes$sample_id,
        receptor = top_clone_celltypes$receptor,
        clone_id_within_sample = top_clone_celltypes$clone_id_within_sample,
        clone_signature = top_clone_celltypes$clone_signature,
        clone_size = top_clone_celltypes$clone_size,
        expansion_category = top_clone_celltypes$expansion_category,
        celltype_label = ifelse(is.na(top_clone_celltypes$celltype_label), "Unknown", top_clone_celltypes$celltype_label)
      ),
      FUN = sum
    )
    names(top_clone_celltypes)[ncol(top_clone_celltypes)] <- "cell_count"
    write.csv(top_clone_celltypes, file.path(out_dir, "top_expanded_clones_with_celltypes.csv"), row.names = FALSE)
  }
}

v_usage <- aggregate(
  rep(1L, nrow(productive_contigs)),
  by = list(
    sample_id = productive_contigs$sample_id,
    receptor = productive_contigs$receptor,
    locus = productive_contigs$locus,
    v_gene = productive_contigs$v_gene
  ),
  FUN = sum
)
names(v_usage)[5] <- "contig_count"
v_usage <- v_usage[v_usage$v_gene != "", , drop = FALSE]
write.csv(v_usage, file.path(out_dir, "v_gene_usage_long.csv"), row.names = FALSE)

j_usage <- aggregate(
  rep(1L, nrow(productive_contigs)),
  by = list(
    sample_id = productive_contigs$sample_id,
    receptor = productive_contigs$receptor,
    locus = productive_contigs$locus,
    j_gene = productive_contigs$j_gene
  ),
  FUN = sum
)
names(j_usage)[5] <- "contig_count"
j_usage <- j_usage[j_usage$j_gene != "", , drop = FALSE]
write.csv(j_usage, file.path(out_dir, "j_gene_usage_long.csv"), row.names = FALSE)

vj_usage <- aggregate(
  rep(1L, nrow(productive_contigs)),
  by = list(
    sample_id = productive_contigs$sample_id,
    receptor = productive_contigs$receptor,
    locus = productive_contigs$locus,
    v_gene = productive_contigs$v_gene,
    j_gene = productive_contigs$j_gene
  ),
  FUN = sum
)
names(vj_usage)[6] <- "contig_count"
vj_usage <- vj_usage[vj_usage$v_gene != "" & vj_usage$j_gene != "", , drop = FALSE]
write.csv(vj_usage, file.path(out_dir, "vj_pair_usage_long.csv"), row.names = FALSE)

cdr3_length_groups <- split(
  productive_contigs$cdr3_aa_length,
  interaction(
    productive_contigs$sample_id,
    productive_contigs$receptor,
    productive_contigs$locus,
    drop = TRUE
  )
)
cdr3_length_summary <- do.call(rbind, lapply(names(cdr3_length_groups), function(key) {
  parts <- strsplit(key, "\\.")[[1]]
  vals <- cdr3_length_groups[[key]]
  data.frame(
    sample_id = parts[1],
    receptor = parts[2],
    locus = parts[3],
    n = sum(!is.na(vals)),
    median_cdr3_aa_length = median(vals, na.rm = TRUE),
    mean_cdr3_aa_length = mean(vals, na.rm = TRUE),
    p90_cdr3_aa_length = unname(quantile(vals, probs = 0.9, na.rm = TRUE)),
    stringsAsFactors = FALSE
  )
}))
names(cdr3_length_summary)[4:7] <- c("n", "median_cdr3_aa_length", "mean_cdr3_aa_length", "p90_cdr3_aa_length")
write.csv(cdr3_length_summary, file.path(out_dir, "cdr3_length_summary.csv"), row.names = FALSE)

top_cdr3 <- aggregate(
  rep(1L, nrow(productive_contigs)),
  by = list(
    sample_id = productive_contigs$sample_id,
    receptor = productive_contigs$receptor,
    locus = productive_contigs$locus,
    cdr3_aa = productive_contigs$cdr3_aa
  ),
  FUN = sum
)
names(top_cdr3)[5] <- "contig_count"
top_cdr3 <- top_cdr3[order(top_cdr3$receptor, -top_cdr3$contig_count, top_cdr3$sample_id), , drop = FALSE]
top_cdr3 <- do.call(rbind, lapply(split(top_cdr3, interaction(top_cdr3$sample_id, top_cdr3$receptor, drop = TRUE)), function(df) {
  head(df, 20)
}))
write.csv(top_cdr3, file.path(out_dir, "top_cdr3_sequences_by_sample.csv"), row.names = FALSE)

plot_sample_order <- sample_ids
sample_receptor_summary$sample_id <- factor(sample_receptor_summary$sample_id, levels = plot_sample_order)
diversity_plot <- ggplot(sample_receptor_summary, aes(x = sample_id, y = shannon, fill = receptor)) +
  geom_col(position = "dodge") +
  labs(x = "Sample", y = "Shannon diversity", title = "Sample-level clonotype diversity") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(out_dir, "sample_clonotype_diversity.png"), diversity_plot, width = 12, height = 5, dpi = 300)
ggsave(file.path(out_dir, "sample_clonotype_diversity.pdf"), diversity_plot, width = 12, height = 5)

paired_plot <- ggplot(sample_receptor_summary, aes(x = sample_id, y = paired_fraction, fill = receptor)) +
  geom_col(position = "dodge") +
  labs(x = "Sample", y = "Paired-chain fraction", title = "Productive paired-chain fraction by sample") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(out_dir, "paired_chain_fraction_by_sample.png"), paired_plot, width = 12, height = 5, dpi = 300)
ggsave(file.path(out_dir, "paired_chain_fraction_by_sample.pdf"), paired_plot, width = 12, height = 5)

expansion_plot <- ggplot(sample_receptor_summary, aes(x = sample_id, y = max_clone_fraction, fill = receptor)) +
  geom_col(position = "dodge") +
  labs(x = "Sample", y = "Largest clone fraction", title = "Largest clone fraction by sample") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(out_dir, "largest_clone_fraction_by_sample.png"), expansion_plot, width = 12, height = 5, dpi = 300)
ggsave(file.path(out_dir, "largest_clone_fraction_by_sample.pdf"), expansion_plot, width = 12, height = 5)

overall_v <- aggregate(
  contig_count ~ receptor + v_gene,
  data = v_usage,
  FUN = sum
)
overall_v <- overall_v[order(overall_v$receptor, -overall_v$contig_count), , drop = FALSE]
overall_v <- do.call(rbind, lapply(split(overall_v, overall_v$receptor), function(df) head(df, 15)))
overall_v_plot <- ggplot(overall_v, aes(x = reorder(v_gene, contig_count), y = contig_count, fill = receptor)) +
  geom_col() +
  facet_wrap(~ receptor, scales = "free_y") +
  coord_flip() +
  labs(x = "V gene", y = "Productive contigs", title = "Top V genes across the 23-sample VDJ set") +
  theme_bw(base_size = 11)
ggsave(file.path(out_dir, "top_v_genes_overall.png"), overall_v_plot, width = 10, height = 7, dpi = 300)
ggsave(file.path(out_dir, "top_v_genes_overall.pdf"), overall_v_plot, width = 10, height = 7)

overall_linked <- if ("rna_metadata_found" %in% names(wide_cell)) sum(wide_cell$rna_metadata_found, na.rm = TRUE) else NA_integer_
overall_cells <- nrow(wide_cell)
overall_link_frac <- if (!is.na(overall_linked) && overall_cells > 0) overall_linked / overall_cells else NA_real_
sample_link_range <- if (exists("link_summary") && nrow(link_summary)) {
  range(link_summary$linkage_fraction, na.rm = TRUE)
} else {
  c(NA_real_, NA_real_)
}
dup_pairs <- sum(file_audit$identical_pair, na.rm = TRUE)
jyh_issue <- file_audit[file_audit$sample_id == "V08", "audit_issue"]

report_lines <- c(
  "# 23例 VDJ 分析小结",
  "",
  paste0("生成时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## 数据源审计",
  paste0("- 读取样本数：", nrow(file_audit), "。"),
  paste0("- BCR/TCR 两目录中有 ", dup_pairs, "/23 对样本 AIRR 文件 MD5 完全一致，说明目录标签不能直接视为独立 BCR 与 TCR 结果。"),
  paste0("- `V08` 的 BCR 目录 AIRR 文件不可读（存在空字节），本次统一使用可读的 TCR 目录文件作为该样本的标准输入。"),
  "- 因此本次分析按 `locus` 字段重新拆分 `BCR`/`TCR`，而不是按目录名拆分。",
  "",
  "## 方法",
  "- 以样本级 AIRR `Dominant_Contigs_AIRR.tsv` 为输入，先审计文件可读性和 BCR/TCR 目录重复情况，再按 `locus` 归并为 BCR/TCR。",
  "- 仅保留 `productive` 且带 `cdr3_aa` 的 contig 进入细胞级克隆分析。",
  paste0("- 用 `sample_id + raw_barcode` 将 VDJ 细胞表与 RNA 元数据关联；本轮 RNA 元数据源为 `", metadata_source, "`。"),
  "- 配对链定义为 BCR 的 `IGH + IGK/IGL` 或 TCR 的 `TRA + TRB`。",
  "",
  "## 主要产出",
  paste0("- 细胞级接口表：`vdj_cell_level_with_metadata.tsv.gz`，可按 `sample_id + raw_barcode` 与 RNA 元数据继续关联。"),
  paste0("- 当前 RNA 元数据来源：`", metadata_source, "`。"),
  paste0("- 有任一 productive 受体的细胞数：", overall_cells, "。"),
  if (!is.na(overall_linked)) paste0(
    "- 其中可直接链接到 RNA 元数据的细胞数：", overall_linked,
    "（总体 ", sprintf("%.1f%%", 100 * overall_link_frac),
    "；样本范围 ", sprintf("%.1f%%", 100 * sample_link_range[1]),
    " - ", sprintf("%.1f%%", 100 * sample_link_range[2]), "）。"
  ) else "- 本轮未找到可用 RNA 元数据文件，因此未做细胞级联表。",
  paste0("- 样本级摘要：`sample_receptor_summary.csv`。"),
  paste0("- 克隆与细胞类型关系：`clone_relation_by_celltype.csv`", if (!is.na(celltype_field)) paste0("（基于字段 `", celltype_field, "`）") else "（本轮无细胞类型字段，未生成）", "。"),
  "- V/J 与 CDR3 结果见 `v_gene_usage_long.csv`、`j_gene_usage_long.csv`、`vj_pair_usage_long.csv`、`top_cdr3_sequences_by_sample.csv`。",
  "",
  "## 已知限制",
  "- 本批 AIRR 文件是混合 VDJ 结果；BCR/TCR 目录多数为重复拷贝，后续若收到官方更正目录，应重新比对。",
  "- 细胞类型关联目前优先使用现有 RNA 元数据中的字段；若后续生成新的重聚类/重注释 metadata，本脚本会优先读取更高层级结果文件。",
  "- 本步骤使用 RNA 元数据连接 VDJ；当前矩阵中的 ADT/WNN 已另行重建，旧供应商 RDS/WNN 标签不用于本次连接。"
)
writeLines(report_lines, con = file.path(out_dir, "VDJ分析报告.md"))

writeLines(capture.output(sessionInfo()), con = file.path(out_dir, "sessionInfo.txt"))
