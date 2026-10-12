#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(scales)
  library(scattermore)
  library(ragg)
  library(grid)
})

set.seed(20260728)

workspace <- Sys.getenv("MPP_WORKSPACE", "inputs/workspace")
package_dir <- file.path(workspace, "工作文件", "论文主图_2026-07-28")
data_dir <- file.path(package_dir, "data")
figure_dir <- file.path(package_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

discovery_root <- file.path(
  workspace,
  "BD_Rhapsody_setup/mac_portable_full/bd_rhapsody_10/downstream",
  "10samples_scrna_mex_v4/results"
)
cross_root <- file.path(
  workspace,
  "BD_Rhapsody_setup/figure_work_package/cross_cohort_results"
)

input_paths <- c(
  discovery_atlas = file.path(data_dir, "discovery_atlas_umap.csv.gz"),
  validation_atlas = file.path(data_dir, "validation_atlas_umap.csv.gz"),
  discovery_myeloid = file.path(data_dir, "discovery_myeloid_umap.csv.gz"),
  validation_myeloid = file.path(data_dir, "validation_myeloid_umap.csv.gz"),
  qc_retention = file.path(data_dir, "qc_retention.csv"),
  marker_expression = file.path(data_dir, "canonical_marker_pseudobulk_expression.csv"),
  pathway_loo_summary = file.path(data_dir, "pathway_loo_summary.csv"),
  pathway_scores = file.path(
    cross_root,
    "04_replicated_mechanism/replicated_pathway_scores_by_sample.csv"
  ),
  pathway_models = file.path(
    cross_root,
    "04_replicated_mechanism/replicated_pathway_score_models.csv"
  ),
  pathway_meta = file.path(
    cross_root,
    "04_replicated_mechanism/replicated_pathway_score_meta_analysis.csv"
  ),
  myeloid_state_composition = file.path(
    cross_root,
    "04_replicated_mechanism/harmonized_macrophage_state_composition.csv"
  ),
  myeloid_state_meta = file.path(
    cross_root,
    "04_replicated_mechanism/harmonized_macrophage_state_meta_analysis.csv"
  ),
  discovery_sensitivity = file.path(
    discovery_root,
    "04_compare/pseudobulk_sensitivity_concordance.csv"
  ),
  validation_concordance = file.path(
    cross_root,
    "03_meta_pathway/validation_adjustment_concordance.csv"
  ),
  composition_meta = file.path(
    cross_root,
    "03_meta_pathway/two_cohort_composition_meta_analysis.csv"
  )
)

missing_inputs <- input_paths[!file.exists(input_paths)]
if (length(missing_inputs)) {
  stop("Missing figure inputs. Run 01_prepare_derived_data.R first:\n",
       paste(missing_inputs, collapse = "\n"))
}

font_family <- "Times New Roman"
set3 <- c(
  "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462",
  "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F"
)
cohort_colors <- c(Discovery = set3[5], Validation = set3[6])
severity_colors <- c(Mild = set3[5], Severe = set3[4])
lineage_colors <- c(
  "B lineage" = set3[3],
  "Macrophage/monocyte" = set3[1],
  "Neutrophil" = set3[7],
  "T lineage" = set3[8],
  "Other" = set3[9]
)
state_colors <- c(
  "Alveolar macrophages" = set3[1],
  "Cycling myeloid" = set3[6],
  "Other myeloid" = set3[9]
)

theme_pub <- theme_classic(base_size = 8.5, base_family = font_family) +
  theme(
    plot.title = element_text(size = 9.5, face = "bold", hjust = 0),
    axis.title = element_text(size = 8.5),
    axis.text = element_text(size = 7.5, color = "black"),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7.5),
    legend.key.height = unit(3.5, "mm"),
    legend.key.width = unit(4.5, "mm"),
    strip.background = element_blank(),
    strip.text = element_text(size = 8, face = "bold", color = "grey20"),
    panel.spacing = unit(3, "mm"),
    plot.margin = margin(4, 5, 4, 5)
  )

tag_theme <- theme(
  plot.tag = element_text(
    family = font_family,
    face = "bold",
    size = 11,
    hjust = 0,
    vjust = 1
  ),
  plot.tag.position = c(0, 1),
  plot.margin = margin(14, 7, 6, 14)
)

save_figure <- function(plot, stem, width, height) {
  pdf_path <- file.path(figure_dir, paste0(stem, ".pdf"))
  png_path <- file.path(figure_dir, paste0(stem, ".png"))
  svg_path <- file.path(figure_dir, paste0(stem, ".svg"))
  tiff_path <- file.path(figure_dir, paste0(stem, ".tiff"))
  replace_pdf_fonts <- function(grob) {
    if (!is.null(grob$gp) && !is.null(grob$gp$fontfamily)) {
      grob$gp$fontfamily <- "Times"
    }
    if (!is.null(grob$children)) {
      grob$children <- lapply(grob$children, replace_pdf_fonts)
    }
    if (!is.null(grob$grobs)) {
      grob$grobs <- lapply(grob$grobs, replace_pdf_fonts)
    }
    grob
  }
  grDevices::pdf(
    file = pdf_path,
    width = width,
    height = height,
    family = "Times",
    useDingbats = FALSE,
    compress = TRUE
  )
  grid::grid.newpage()
  grid::grid.draw(replace_pdf_fonts(ggplot2::ggplotGrob(plot)))
  grDevices::dev.off()
  ggsave(
    png_path,
    plot = plot,
    device = ragg::agg_png,
    width = width,
    height = height,
    units = "in",
    dpi = 300,
    bg = "white"
  )
  ggsave(
    svg_path,
    plot = plot,
    device = svglite::svglite,
    width = width,
    height = height,
    units = "in",
    bg = "white"
  )
  ggsave(
    tiff_path,
    plot = plot,
    device = ragg::agg_tiff,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    compression = "lzw",
    bg = "white"
  )
  invisible(c(pdf = pdf_path, png = png_path, svg = svg_path, tiff = tiff_path))
}

pathway_labels <- c(
  "GO:0002250" = "Adaptive immune response",
  "GO:0009615" = "Response to virus",
  "GO:0031295" = "T-cell costimulation",
  "GO:0033209" = "TNF-mediated signaling",
  "GO:0050862" = "Positive regulation of TCR signaling",
  "GO:0050871" = "Positive regulation of B-cell activation",
  "GO:0071222" = "Cellular response to LPS"
)
pathway_order <- unname(pathway_labels)

display_label_audit <- c(
  names(lineage_colors), names(state_colors), pathway_order,
  "Discovery cohort", "Validation cohort", "Mild", "Severe"
)
if (any(grepl("_|GO:", display_label_audit))) {
  stop("Display label audit failed: underscore or raw GO identifier remains")
}

# Figure 1: study design, QC, harmonized atlas, and canonical markers.
design_boxes <- data.table(
  xmin = c(0.2, 5.25, 1.3, 5.3),
  xmax = c(4.75, 9.8, 4.7, 9.65),
  ymin = c(1.6, 1.6, 0.1, 0.1),
  ymax = c(2.85, 2.85, 1.0, 1.0),
  fill = c("Discovery", "Validation", "Harmonized", "Mechanism"),
  label = c(
    "Discovery cohort\n10 patients | 5 mild, 5 severe\n115,594 cells",
    "Validation cohort\n23 patients | 10 mild, 13 severe\n209,299 cells",
    "Matched quality control and\nlineage annotation",
    "Patient-level pathway analysis and\ncross-cohort synthesis"
  )
)
p1a <- ggplot() +
  geom_rect(
    data = design_boxes,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "white",
    color = "grey35",
    linewidth = 0.7
  ) +
  geom_text(
    data = design_boxes,
    aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = label),
    family = font_family,
    size = 2.65,
    lineheight = 1.1
  ) +
  annotate(
    "segment",
    x = c(2.5, 7.5, 4.7),
    xend = c(2.8, 3.7, 5.3),
    y = c(1.58, 1.58, 0.55),
    yend = c(1.02, 1.02, 0.55),
    linewidth = 0.5,
    color = "grey25",
    arrow = arrow(type = "closed", length = unit(2.3, "mm"))
  ) +
  coord_cartesian(xlim = c(0, 10), ylim = c(0, 3), clip = "off") +
  labs(title = "Study design and cross-cohort analysis") +
  theme_void(base_family = font_family) +
  theme(
    plot.title = element_text(size = 9.5, face = "bold", hjust = 0),
    legend.position = "none",
    plot.margin = margin(3, 4, 2, 4)
  )

qc <- fread(input_paths[["qc_retention"]])
qc[, cohort := factor(cohort, levels = c("Discovery", "Validation"))]
qc[, severity_label := factor(
  fifelse(severity == "mild", "Mild", "Severe"),
  levels = c("Mild", "Severe")
)]
qc[, group := factor(
  paste(cohort, severity_label, sep = "\n"),
  levels = c(
    "Discovery\nMild", "Discovery\nSevere",
    "Validation\nMild", "Validation\nSevere"
  )
)]
p1b <- ggplot(qc, aes(x = severity_label, y = retention_rate, color = severity_label)) +
  geom_boxplot(
    aes(group = severity_label),
    width = 0.55,
    outlier.shape = NA,
    fill = "white",
    linewidth = 0.45
  ) +
  geom_jitter(width = 0.12, height = 0, size = 1.45, alpha = 0.9) +
  scale_color_manual(values = severity_colors) +
  facet_grid(. ~ cohort) +
  scale_y_continuous(
    limits = c(0, 1.02),
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = "Cell retention after QC",
    x = NULL,
    y = "Retained cells (%)",
    color = NULL
  ) +
  theme_pub +
  theme(
    axis.text.x = element_text(size = 7),
    strip.text = element_text(size = 7.5),
    legend.position = "none"
  )

make_atlas_umap <- function(path, title_text, expected_cells) {
  atlas <- fread(path)
  if (nrow(atlas) != expected_cells) {
    stop("Unexpected atlas cell count for ", title_text)
  }
  atlas[, lineage := factor(lineage, levels = names(lineage_colors))]
  atlas[, draw_order := lineage != "Other"]
  setorder(atlas, draw_order)
  label_data <- atlas[
    lineage != "Other",
    .(umap_1 = median(umap_1), umap_2 = median(umap_2)),
    by = lineage
  ]
  label_data[, display_label := fcase(
    lineage == "Macrophage/monocyte", "Macrophages and\nmonocytes",
    default = as.character(lineage)
  )]
  ggplot(atlas, aes(x = umap_1, y = umap_2, color = lineage)) +
    scattermore::geom_scattermore(
      pointsize = 1.2,
      pixels = c(1200, 1200),
      alpha = 0.75
    ) +
    geom_text_repel(
      data = label_data,
      aes(label = display_label),
      family = font_family,
      size = 2.3,
      color = "black",
      box.padding = unit(0.45, "lines"),
      point.padding = unit(1.6, "mm"),
      force = 1.25,
      max.overlaps = Inf,
      min.segment.length = 0,
      segment.color = alpha("grey35", 0.65),
      segment.size = 0.3,
      seed = 20260728,
      show.legend = FALSE
    ) +
    scale_color_manual(values = lineage_colors, drop = FALSE) +
    coord_equal() +
    labs(title = title_text, x = "UMAP 1", y = "UMAP 2", color = NULL) +
    theme_pub +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "none"
    )
}

p1c <- make_atlas_umap(
  input_paths[["discovery_atlas"]],
  "Discovery atlas",
  115594
)
p1d <- make_atlas_umap(
  input_paths[["validation_atlas"]],
  "Validation atlas",
  209299
)

marker_expression <- fread(input_paths[["marker_expression"]])
marker_gene_order <- c(
  "MS4A1", "CD79A", "CD74",
  "LST1", "FCER1G", "CTSS",
  "FCGR3B", "CSF3R", "S100A8",
  "CD3D", "CD3E", "TRBC1"
)
marker_expression[, gene := factor(gene, levels = rev(marker_gene_order))]
marker_expression[, lineage := factor(
  lineage,
  levels = c("B lineage", "Macrophage/monocyte", "Neutrophil", "T lineage")
)]
marker_expression[, cohort := factor(cohort, levels = c("Discovery", "Validation"))]
p1e <- ggplot(marker_expression, aes(x = lineage, y = gene, fill = expression_z)) +
  geom_tile(color = "white", linewidth = 0.35) +
  scale_fill_gradient2(
    low = set3[5],
    mid = "white",
    high = set3[4],
    midpoint = 0,
    limits = c(-2.3, 2.3),
    oob = squish,
    name = "Scaled\nexpression"
  ) +
  facet_grid(. ~ cohort) +
  scale_x_discrete(labels = c(
    "B lineage" = "B lineage",
    "Macrophage/monocyte" = "Macrophage and\nmonocyte",
    "Neutrophil" = "Neutrophil",
    "T lineage" = "T lineage"
  )) +
  labs(
    title = "Canonical lineage markers",
    x = NULL,
    y = NULL
  ) +
  theme_pub +
  theme(
    axis.text.x = element_text(size = 6.3, lineheight = 0.9),
    axis.text.y = element_text(face = "italic"),
    panel.grid = element_blank(),
    legend.position = "right"
  )

figure_1 <- p1a /
  (p1b | p1c | p1d) /
  p1e +
  plot_layout(heights = c(1.25, 2.35, 2.15)) +
  plot_annotation(tag_levels = "A", theme = tag_theme)
save_figure(figure_1, "Figure_1_study_design_QC_atlas", 6.9, 8.15)

# Figure 2: seven replicated pathways and cohort/meta forest plot.
pathway_scores <- fread(input_paths[["pathway_scores"]])
pathway_scores[, pathway := unname(pathway_labels[go_id])]
if (anyNA(pathway_scores$pathway)) {
  stop("Unmapped pathway label in pathway score table")
}
pathway_scores[, pathway := factor(pathway, levels = rev(pathway_order))]
pathway_scores[, cohort_label := factor(
  fcase(cohort == "discovery_10", "Discovery (n = 9)",
        cohort == "validation_23", "Validation (n = 23)"),
  levels = c("Discovery (n = 9)", "Validation (n = 23)")
)]
pathway_scores[, severity_label := factor(
  fifelse(severity == "mild", "Mild", "Severe"),
  levels = c("Mild", "Severe")
)]
sample_order <- unique(pathway_scores[, .(cohort_label, sample_id, severity_label)])
setorder(sample_order, cohort_label, severity_label, sample_id)
sample_order[, patient_index := seq_len(.N), by = cohort_label]
pathway_scores <- merge(
  pathway_scores,
  sample_order,
  by = c("cohort_label", "sample_id", "severity_label"),
  sort = FALSE
)
severity_strip <- unique(pathway_scores[, .(cohort_label, patient_index, severity_label)])
severity_strip[, pathway_position := 8]
pathway_scores[, pathway_position := as.numeric(pathway)]
p2a <- ggplot(
  pathway_scores,
  aes(x = patient_index, y = pathway_position, fill = standardized_score)
) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_point(
    data = severity_strip,
    aes(x = patient_index, y = pathway_position, color = severity_label),
    inherit.aes = FALSE,
    shape = 15,
    size = 1.8
  ) +
  facet_grid(. ~ cohort_label, scales = "free_x", space = "free_x") +
  scale_y_continuous(
    breaks = seq_along(rev(pathway_order)),
    labels = rev(pathway_order),
    limits = c(0.5, 8.35),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  scale_fill_gradient2(
    low = set3[5],
    mid = "white",
    high = set3[4],
    midpoint = 0,
    limits = c(-2.5, 2.5),
    oob = squish,
    name = "Standardized\npathway score"
  ) +
  scale_color_manual(values = severity_colors, name = NULL) +
  labs(
    title = "Patient-level pathway scores",
    x = "Patients",
    y = NULL
  ) +
  theme_pub +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.box = "vertical"
  )

pathway_meta <- fread(input_paths[["pathway_meta"]])
pathway_forest <- rbindlist(list(
  pathway_meta[, .(
    go_id,
    estimate = estimate_discovery,
    se = standard_error_discovery,
    estimate_type = "Discovery"
  )],
  pathway_meta[, .(
    go_id,
    estimate = estimate_validation,
    se = standard_error_validation,
    estimate_type = "Validation"
  )],
  pathway_meta[, .(
    go_id,
    estimate = fixed_effect,
    se = fixed_se,
    estimate_type = "Fixed-effect meta"
  )],
  pathway_meta[, .(
    go_id,
    estimate = random_effect,
    se = random_se,
    estimate_type = "Random-effects meta"
  )]
))
pathway_forest[, pathway := unname(pathway_labels[go_id])]
pathway_forest[, pathway_base := match(pathway, rev(pathway_order))]
forest_offsets <- c(
  Discovery = 0.27,
  Validation = 0.09,
  `Fixed-effect meta` = -0.09,
  `Random-effects meta` = -0.27
)
pathway_forest[, y := pathway_base + unname(forest_offsets[estimate_type])]
pathway_forest[, `:=`(
  ci_low = estimate - 1.96 * se,
  ci_high = estimate + 1.96 * se,
  estimate_type = factor(
    estimate_type,
    levels = c("Discovery", "Validation", "Fixed-effect meta", "Random-effects meta")
  )
)]
forest_colors <- c(
  Discovery = cohort_colors[["Discovery"]],
  Validation = cohort_colors[["Validation"]],
  `Fixed-effect meta` = "black",
  `Random-effects meta` = "grey40"
)
forest_shapes <- c(
  Discovery = 16,
  Validation = 17,
  `Fixed-effect meta` = 18,
  `Random-effects meta` = 23
)
p2b <- ggplot(pathway_forest, aes(x = estimate, y = y, color = estimate_type)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.45) +
  geom_errorbar(
    aes(xmin = ci_low, xmax = ci_high),
    orientation = "y",
    width = 0,
    linewidth = 0.55
  ) +
  geom_point(aes(shape = estimate_type), size = 2.25, stroke = 0.8) +
  scale_color_manual(values = forest_colors, name = NULL) +
  scale_shape_manual(values = forest_shapes, name = NULL) +
  scale_y_continuous(
    breaks = seq_along(rev(pathway_order)),
    labels = rev(pathway_order),
    expand = expansion(add = c(0.55, 0.55))
  ) +
  labs(
    title = "Cohort and meta-analysis effects",
    x = "Severe - mild standardized score (95% CI)",
    y = NULL
  ) +
  theme_pub +
  theme(
    panel.grid.major.y = element_line(color = "grey96", linewidth = 0.3),
    legend.position = "bottom",
    legend.box.spacing = unit(1, "mm")
  ) +
  guides(color = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))

figure_2 <- p2a / p2b +
  plot_layout(heights = c(2.55, 3.7)) +
  plot_annotation(tag_levels = "A", theme = tag_theme)
save_figure(figure_2, "Figure_2_replicated_pathways", 6.9, 7.25)

# Figure 3: harmonized myeloid states and exploratory composition analysis.
make_state_umap <- function(path, title_text, expected_cells) {
  state_data <- fread(path)
  if (nrow(state_data) != expected_cells) {
    stop("Unexpected myeloid-state cell count for ", title_text)
  }
  state_data[, state := factor(
    state,
    levels = c("Other myeloid", "Alveolar macrophages", "Cycling myeloid")
  )]
  background_data <- state_data[state == "Other myeloid"]
  highlight_data <- state_data[state != "Other myeloid"]
  legend_data <- data.table(
    umap_1 = NA_real_,
    umap_2 = NA_real_,
    state = factor(
      c("Other myeloid", "Alveolar macrophages", "Cycling myeloid"),
      levels = levels(state_data$state)
    )
  )
  ggplot() +
    scattermore::geom_scattermore(
      data = background_data,
      aes(x = umap_1, y = umap_2, color = state),
      pointsize = 0.9,
      pixels = c(1200, 1200),
      alpha = 0.38
    ) +
    scattermore::geom_scattermore(
      data = highlight_data,
      aes(x = umap_1, y = umap_2, color = state),
      pointsize = 1.8,
      pixels = c(1200, 1200),
      alpha = 0.92
    ) +
    geom_point(
      data = legend_data,
      aes(x = umap_1, y = umap_2, color = state),
      shape = 16,
      size = 2.4,
      alpha = 0,
      inherit.aes = FALSE,
      show.legend = TRUE
    ) +
    scale_color_manual(values = state_colors, drop = FALSE) +
    coord_equal() +
    labs(title = title_text, x = "UMAP 1", y = "UMAP 2", color = NULL) +
    theme_pub +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "bottom"
    ) +
    guides(color = guide_legend(
      override.aes = list(shape = 16, size = 2.4, alpha = 1),
      nrow = 1
    ))
}

p3a <- make_state_umap(
  input_paths[["discovery_myeloid"]],
  "Discovery cohort",
  23228
)
p3b <- make_state_umap(
  input_paths[["validation_myeloid"]],
  "Validation cohort",
  118518
) + theme(legend.position = "none")

state_composition <- fread(input_paths[["myeloid_state_composition"]])
state_composition[, cohort_label := fcase(
  cohort == "discovery_10", "Discovery",
  cohort == "validation_23", "Validation"
)]
state_composition[, severity_label := factor(
  fifelse(severity == "mild", "Mild", "Severe"),
  levels = c("Mild", "Severe")
)]
state_composition[, group := factor(
  paste(cohort_label, severity_label, sep = "\n"),
  levels = c(
    "Discovery\nMild", "Discovery\nSevere",
    "Validation\nMild", "Validation\nSevere"
  )
)]

make_state_fraction_plot <- function(state_id, title_text, y_limit) {
  plot_data <- state_composition[harmonized_state == state_id]
  ggplot(plot_data, aes(x = group, y = fraction, color = severity_label)) +
    geom_boxplot(
      aes(group = group),
      width = 0.55,
      outlier.shape = NA,
      fill = "white",
      linewidth = 0.45
    ) +
    geom_jitter(width = 0.12, height = 0, size = 1.45, alpha = 0.9) +
    scale_color_manual(values = severity_colors) +
    scale_y_continuous(
      limits = c(0, y_limit),
      labels = percent_format(accuracy = 1),
      expand = expansion(mult = c(0, 0.03))
    ) +
    labs(title = title_text, x = NULL, y = "Fraction of myeloid cells", color = NULL) +
    theme_pub +
    theme(
      axis.text.x = element_text(size = 7, lineheight = 0.95),
      legend.position = "none"
    )
}

p3c <- make_state_fraction_plot(
  "Alveolar_macrophage",
  "Alveolar macrophages",
  0.9
)
p3d <- make_state_fraction_plot(
  "Cycling_myeloid",
  "Cycling myeloid",
  0.25
)

state_meta <- fread(input_paths[["myeloid_state_meta"]])
state_forest <- rbindlist(list(
  state_meta[, .(
    harmonized_state,
    estimate = estimate_log_odds_discovery,
    se = standard_error_discovery,
    estimate_type = "Discovery"
  )],
  state_meta[, .(
    harmonized_state,
    estimate = estimate_log_odds_validation,
    se = standard_error_validation,
    estimate_type = "Validation"
  )],
  state_meta[, .(
    harmonized_state,
    estimate = fixed_effect,
    se = fixed_se,
    estimate_type = "Fixed-effect meta"
  )]
))
state_forest[, state := fcase(
  harmonized_state == "Alveolar_macrophage", "Alveolar macrophages",
  harmonized_state == "Cycling_myeloid", "Cycling myeloid"
)]
state_forest[, state_base := fcase(
  state == "Alveolar macrophages", 2,
  state == "Cycling myeloid", 1
)]
state_offsets <- c(Discovery = 0.2, Validation = 0, `Fixed-effect meta` = -0.2)
state_forest[, y := state_base + unname(state_offsets[estimate_type])]
state_forest[, `:=`(
  ci_low = estimate - 1.96 * se,
  ci_high = estimate + 1.96 * se,
  estimate_type = factor(
    estimate_type,
    levels = c("Discovery", "Validation", "Fixed-effect meta")
  )
)]
p3e <- ggplot(state_forest, aes(x = estimate, y = y, color = estimate_type)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.45) +
  geom_errorbar(
    aes(xmin = ci_low, xmax = ci_high),
    orientation = "y",
    width = 0,
    linewidth = 0.6
  ) +
  geom_point(aes(shape = estimate_type), size = 2.35) +
  scale_color_manual(
    values = c(
      Discovery = cohort_colors[["Discovery"]],
      Validation = cohort_colors[["Validation"]],
      `Fixed-effect meta` = "black"
    ),
    name = NULL
  ) +
  scale_shape_manual(
    values = c(Discovery = 16, Validation = 17, `Fixed-effect meta` = 18),
    name = NULL
  ) +
  scale_y_continuous(
    breaks = c(1, 2),
    labels = c("Cycling myeloid", "Alveolar macrophages"),
    expand = expansion(add = c(0.5, 0.5))
  ) +
  labs(
    title = "State-composition effects",
    x = "Severe - mild log odds (95% CI)",
    y = NULL
  ) +
  theme_pub +
  theme(
    panel.grid.major.y = element_line(color = "grey96", linewidth = 0.3),
    legend.position = "bottom"
  ) +
  guides(color = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))

figure_3 <- (p3a | p3b) /
  (p3c | p3d) /
  p3e +
  plot_layout(heights = c(2.5, 2.25, 1.8)) +
  plot_annotation(tag_levels = "A", theme = tag_theme)
save_figure(figure_3, "Figure_3_macrophage_states", 6.9, 8.05)

# Figure 4: covariate, LOO, QC-outlier, and composition robustness.
pathway_models <- fread(input_paths[["pathway_models"]])
validation_effects <- pathway_models[
  cohort == "validation_23" & coefficient == "severitysevere"
]
validation_effects[, model_label := factor(
  fcase(
    model_id == "severity_only", "Primary",
    model_id == "covariate_adjusted", "Covariate adjusted"
  ),
  levels = c("Primary", "Covariate adjusted")
)]
validation_effects[, pathway := factor(
  unname(pathway_labels[go_id]),
  levels = rev(pathway_order)
)]
validation_effects[, `:=`(
  ci_low = estimate - 1.96 * standard_error,
  ci_high = estimate + 1.96 * standard_error
)]
p4a <- ggplot(
  validation_effects,
  aes(x = estimate, y = pathway, color = model_label, group = pathway)
) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.45) +
  geom_line(color = "grey75", linewidth = 0.45) +
  geom_errorbar(
    aes(xmin = ci_low, xmax = ci_high),
    position = position_dodge(width = 0.45),
    orientation = "y",
    width = 0,
    linewidth = 0.55
  ) +
  geom_point(
    aes(shape = model_label),
    position = position_dodge(width = 0.45),
    size = 2.1
  ) +
  scale_color_manual(
    values = c(Primary = set3[5], `Covariate adjusted` = set3[4]),
    name = NULL
  ) +
  scale_shape_manual(values = c(Primary = 16, `Covariate adjusted` = 17), name = NULL) +
  labs(
    title = "Covariate adjustment",
    x = "Severe - mild effect (95% CI)",
    y = NULL
  ) +
  theme_pub +
  theme(
    panel.grid.major.y = element_line(color = "grey96", linewidth = 0.3),
    legend.position = "bottom"
  ) +
  guides(color = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))

discovery_sensitivity <- fread(input_paths[["discovery_sensitivity"]])
discovery_sensitivity <- discovery_sensitivity[
  cell_type_broad %chin% c("B_Plasma", "Macrophage_Monocyte", "Neutrophil", "T_NK") &
    scenario_id == "sensitivity_exclude_all_qc_outliers"
]
discovery_sensitivity[, lineage := fcase(
  cell_type_broad == "B_Plasma", "B lineage",
  cell_type_broad == "Macrophage_Monocyte", "Macrophage/monocyte",
  cell_type_broad == "Neutrophil", "Neutrophil",
  cell_type_broad == "T_NK", "T lineage"
)]
discovery_sensitivity_long <- melt(
  discovery_sensitivity,
  id.vars = "lineage",
  measure.vars = c("logFC_spearman", "direction_agreement"),
  variable.name = "metric",
  value.name = "value"
)
discovery_sensitivity_long[, metric := factor(
  fcase(
    metric == "logFC_spearman", "Effect correlation",
    metric == "direction_agreement", "Sign agreement"
  ),
  levels = c("Effect correlation", "Sign agreement")
)]
p4b <- ggplot(
  discovery_sensitivity_long,
  aes(x = value, y = lineage, color = metric, shape = metric)
) +
  geom_vline(xintercept = 0.75, linetype = "dotted", color = "grey65", linewidth = 0.4) +
  geom_point(position = position_dodge(width = 0.55), size = 2.15) +
  scale_color_manual(
    values = c("Effect correlation" = set3[5], "Sign agreement" = set3[1]),
    name = NULL
  ) +
  scale_shape_manual(values = c("Effect correlation" = 16, "Sign agreement" = 17), name = NULL) +
  scale_x_continuous(
    limits = c(0.6, 1.01),
    breaks = seq(0.6, 1, 0.1),
    labels = number_format(accuracy = 0.1)
  ) +
  labs(
    title = "QC sensitivity",
    x = "Concordance",
    y = NULL
  ) +
  scale_y_discrete(labels = c("Macrophage/monocyte" = "Macrophages and monocytes")) +
  theme_pub +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE), shape = "none")

loo_summary <- fread(input_paths[["pathway_loo_summary"]])
loo_summary[, pathway := factor(
  unname(pathway_labels[go_id]),
  levels = rev(pathway_order)
)]
loo_summary[, model_label := factor(
  fcase(
    model_key == "discovery_primary", "Discovery primary",
    model_key == "validation_primary", "Validation primary",
    model_key == "validation_adjusted", "Validation adjusted"
  ),
  levels = c("Discovery primary", "Validation primary", "Validation adjusted")
)]
p4c <- ggplot(loo_summary, aes(y = pathway, color = model_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.45) +
  geom_errorbar(
    aes(xmin = min_loo_effect, xmax = max_loo_effect, x = full_effect),
    position = position_dodge(width = 0.58),
    orientation = "y",
    width = 0,
    linewidth = 0.7
  ) +
  geom_point(
    aes(x = full_effect, shape = model_label),
    position = position_dodge(width = 0.58),
    size = 2.1
  ) +
  scale_color_manual(
    values = c(
      `Discovery primary` = set3[5],
      `Validation primary` = set3[6],
      `Validation adjusted` = set3[4]
    ),
    name = NULL
  ) +
  scale_shape_manual(values = c(16, 17, 15), name = NULL) +
  labs(
    title = "Leave-one-patient-out analysis",
    x = "Severe - mild pathway score difference",
    y = NULL
  ) +
  theme_pub +
  theme(
    panel.grid.major.y = element_line(color = "grey96", linewidth = 0.3),
    legend.position = "bottom"
  ) +
  guides(color = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))

composition_meta <- fread(input_paths[["composition_meta"]])
composition_meta[, lineage := factor(
  fcase(
    harmonized_lineage == "B_lineage", "B lineage",
    harmonized_lineage == "Macrophage_Monocyte", "Macrophages and monocytes",
    harmonized_lineage == "Neutrophil", "Neutrophil",
    harmonized_lineage == "T_lineage", "T lineage"
  ),
  levels = rev(c("B lineage", "Macrophages and monocytes", "Neutrophil", "T lineage"))
)]
composition_meta[, `:=`(
  ci_low = fixed_effect - 1.96 * fixed_se,
  ci_high = fixed_effect + 1.96 * fixed_se
)]
p4d <- ggplot(composition_meta, aes(x = fixed_effect, y = lineage)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.45) +
  geom_errorbar(
    aes(xmin = ci_low, xmax = ci_high),
    orientation = "y",
    width = 0,
    linewidth = 0.65
  ) +
  geom_point(shape = 18, size = 2.5, color = "black") +
  labs(
    title = "Broad-lineage composition",
    x = "Severe - mild cell fraction (95% CI)",
    y = NULL
  ) +
  theme_pub +
  theme(panel.grid.major.y = element_line(color = "grey96", linewidth = 0.3))

figure_4 <- (p4a | p4b) /
  p4c /
  p4d +
  plot_layout(heights = c(2.5, 2.7, 1.75)) +
  plot_annotation(tag_levels = "A", theme = tag_theme)
save_figure(figure_4, "Figure_4_robustness", 6.9, 8.2)

figure_checks <- data.table(
  check = c(
    "figure_1_pdf_exists",
    "figure_2_pdf_exists",
    "figure_3_pdf_exists",
    "figure_4_pdf_exists",
    "figure_1_png_exists",
    "figure_2_png_exists",
    "figure_3_png_exists",
    "figure_4_png_exists",
    "figure_1_svg_exists",
    "figure_2_svg_exists",
    "figure_3_svg_exists",
    "figure_4_svg_exists",
    "figure_1_tiff_exists",
    "figure_2_tiff_exists",
    "figure_3_tiff_exists",
    "figure_4_tiff_exists",
    "pathway_score_rows",
    "pathway_meta_rows",
    "pathway_forest_rows",
    "state_composition_rows",
    "state_meta_rows",
    "loo_summary_rows",
    "loo_all_negative",
    "display_label_underscores_or_GO_ids"
  ),
  observed = c(
    file.exists(file.path(figure_dir, "Figure_1_study_design_QC_atlas.pdf")),
    file.exists(file.path(figure_dir, "Figure_2_replicated_pathways.pdf")),
    file.exists(file.path(figure_dir, "Figure_3_macrophage_states.pdf")),
    file.exists(file.path(figure_dir, "Figure_4_robustness.pdf")),
    file.exists(file.path(figure_dir, "Figure_1_study_design_QC_atlas.png")),
    file.exists(file.path(figure_dir, "Figure_2_replicated_pathways.png")),
    file.exists(file.path(figure_dir, "Figure_3_macrophage_states.png")),
    file.exists(file.path(figure_dir, "Figure_4_robustness.png")),
    file.exists(file.path(figure_dir, "Figure_1_study_design_QC_atlas.svg")),
    file.exists(file.path(figure_dir, "Figure_2_replicated_pathways.svg")),
    file.exists(file.path(figure_dir, "Figure_3_macrophage_states.svg")),
    file.exists(file.path(figure_dir, "Figure_4_robustness.svg")),
    file.exists(file.path(figure_dir, "Figure_1_study_design_QC_atlas.tiff")),
    file.exists(file.path(figure_dir, "Figure_2_replicated_pathways.tiff")),
    file.exists(file.path(figure_dir, "Figure_3_macrophage_states.tiff")),
    file.exists(file.path(figure_dir, "Figure_4_robustness.tiff")),
    nrow(pathway_scores),
    nrow(pathway_meta),
    nrow(pathway_forest),
    nrow(state_composition),
    nrow(state_meta),
    nrow(loo_summary),
    sum(loo_summary$all_negative),
    sum(grepl("_|GO:", display_label_audit))
  ),
  expected = c(
    TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
    TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
    224, 7, 28, 64, 2, 21, 21, 0
  )
)
figure_checks[, status := fifelse(as.character(observed) == as.character(expected), "PASS", "FAIL")]
fwrite(figure_checks, file.path(package_dir, "figure_data_checks.csv"))

if (any(figure_checks$status == "FAIL")) {
  stop("Figure-data validation failed; inspect figure_data_checks.csv")
}

writeLines(
  capture.output(sessionInfo()),
  file.path(package_dir, "sessionInfo.txt")
)

message("All four main figures were rendered and passed data checks")
