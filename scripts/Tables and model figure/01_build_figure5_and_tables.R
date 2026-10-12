#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(grid)
  library(ragg)
  library(svglite)
})

workspace <- Sys.getenv("MPP_WORKSPACE", "inputs/workspace")
source_dir <- file.path(workspace, "工作文件")
output_dir <- file.path(source_dir, "文章1", "中文图表整合审阅稿 2026-08-21")
figure_dir <- file.path(output_dir, "figures")
table_dir <- file.path(output_dir, "tables")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

axis_file <- file.path(source_dir, "文章1_胞葬HIF1A轴结果_2026-08-19.csv")
cross_file <- file.path(source_dir, "两篇文章_关键节点跨组学对照_2026-08-19.csv")
protein_file <- file.path(source_dir, "两篇文章_蛋白组关键结果_2026-08-19.csv")
vdj_file <- file.path(source_dir, "文章1_VDJ结果_2026-08-19.csv")
cellchat_file <- file.path(
  source_dir,
  "23例单细胞重新分析", "results", "06_cellchat_by_sample",
  "pathway_severity_comparison.csv"
)

required_files <- c(axis_file, cross_file, protein_file, vdj_file, cellchat_file)
if (any(!file.exists(required_files))) {
  stop("Missing source files:\n", paste(required_files[!file.exists(required_files)], collapse = "\n"))
}

axis_data <- fread(axis_file, encoding = "UTF-8", check.names = FALSE)
cross_data <- fread(cross_file, encoding = "UTF-8", check.names = FALSE)
protein_data <- fread(protein_file, encoding = "UTF-8", check.names = FALSE)
vdj_data <- fread(vdj_file, encoding = "UTF-8", check.names = FALSE)
cellchat_data <- fread(cellchat_file, encoding = "UTF-8", check.names = FALSE)

font_family <- "Times New Roman"
set3 <- c(
  "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462",
  "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F"
)
theme_pub <- theme_classic(base_size = 8.2, base_family = font_family) +
  theme(
    plot.title = element_text(size = 9.4, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 7.1, color = "grey25", hjust = 0),
    axis.title = element_text(size = 8.0),
    axis.text = element_text(size = 7.0, color = "black"),
    legend.title = element_text(size = 7.2),
    legend.text = element_text(size = 6.8),
    plot.margin = margin(6, 7, 5, 7)
  )

panel_tag_theme <- theme(
  plot.tag = element_text(family = font_family, face = "bold", size = 11),
  plot.tag.position = c(0, 1),
  plot.margin = margin(14, 7, 6, 14)
)

save_figure <- function(plot, stem, width = 6.9, height = 9.4) {
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
    pdf_path,
    width = width,
    height = height,
    family = "Times",
    useDingbats = FALSE,
    compress = TRUE
  )
  grid::grid.newpage()
  grid::grid.draw(replace_pdf_fonts(ggplot2::ggplotGrob(plot)))
  grDevices::dev.off()
  ggsave(png_path, plot, device = ragg::agg_png, width = width, height = height,
        units = "in", dpi = 300, bg = "white")
  ggsave(svg_path, plot, device = svglite::svglite, width = width, height = height,
        units = "in", bg = "white")
  ggsave(tiff_path, plot, device = ragg::agg_tiff, width = width, height = height,
        units = "in", dpi = 600, compression = "lzw", bg = "white")
}

# Panel A: discovery-cohort localization of the three efferocytosis-related scores.
score_order <- c("Efferocytosis", "PHAGOCYTOSIS", "PHAGOCYTOSIS_ENGULFMENT")
cell_order <- c(
  "alveolar Macro", "Mono-derived Macro", "Mono-derived alveolar Macro", "Monocytes",
  "mDC", "Neutrophil", "Eosinophil", "NK", "B", "CD4 T", "CD8 T",
  "Proliferating T", "pDC", "Plasma"
)
cell_labels <- c(
  "alveolar Macro" = "Alveolar-like macrophages",
  "Mono-derived Macro" = "Monocyte-derived macrophages",
  "Mono-derived alveolar Macro" = "Monocyte-derived alveolar-like macrophages",
  "Monocytes" = "Monocytes",
  "mDC" = "mDC", "Neutrophil" = "Neutrophils", "Eosinophil" = "Eosinophils",
  "NK" = "NK cells", "B" = "B cells", "CD4 T" = "CD4 T cells",
  "CD8 T" = "CD8 T cells", "Proliferating T" = "Cycling T cells",
  "pDC" = "pDC", "Plasma" = "Plasma cells"
)
score_labels <- c(
  "Efferocytosis" = "Efferocytosis",
  "PHAGOCYTOSIS" = "Phagocytosis",
  "PHAGOCYTOSIS_ENGULFMENT" = "Engulfment"
)

score_data <- axis_data[
  layer == "efferocytosis_score" & score_name %in% score_order,
  .(cell_type, score_name, score)
]
stopifnot(nrow(score_data) == length(cell_order) * length(score_order))
score_data[, cell_type := factor(cell_type, levels = rev(cell_order))]
score_data[, score_name := factor(score_name, levels = score_order)]

p5a <- ggplot(score_data, aes(x = score_name, y = cell_type, color = score)) +
  geom_hline(yintercept = seq(1.5, length(cell_order) - 0.5, by = 1),
             color = "grey97", linewidth = 0.2) +
  geom_point(aes(size = abs(score)), shape = 16) +
  scale_x_discrete(labels = score_labels, position = "top") +
  scale_y_discrete(labels = cell_labels) +
  scale_color_gradient2(
    low = set3[5], mid = "white", high = set3[4], midpoint = 0,
    limits = c(-0.10, 0.22), oob = scales::squish, name = "Mean score"
  ) +
  scale_size(range = c(1.2, 4.8), limits = c(0, 0.22), guide = "none") +
  labs(
    title = "Discovery-cohort localization",
    subtitle = "10-patient discovery cohort",
    x = NULL, y = NULL
  ) +
  theme_pub +
  theme(
    axis.text.x.top = element_text(angle = 28, hjust = 0, vjust = 0),
    axis.text.y = element_text(size = 6.4),
    axis.line = element_blank(), axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.key.width = unit(19, "mm")
  )

# Panel B: candidate-node logFC lollipop plots. The source-defined discovery
# contrast is shown without matrix cells or decorative module boxes.
candidate_gene_order_plot <- c(
  "HIF1A", "THBS1", "CD36", "MERTK", "AXL", "CD47", "SIRPA",
  "C1QA", "C1QB", "C1QC", "APOE", "TREM2", "CTSD", "LAMP1", "NPC2", "LIPA",
  "MAPK1", "MAPK3", "TGFB1", "SMAD2", "SMAD3", "SMAD4"
)
candidate_cell_order_plot <- c("Mono-derived Macro", "Mono-derived alveolar Macro", "Monocytes")
candidate_module_plot <- c(
  HIF1A = "Hypoxia and entry",
  THBS1 = "Recognition and contact", CD36 = "Recognition and contact",
  MERTK = "Recognition and contact", AXL = "Recognition and contact",
  CD47 = "Engulfment brake", SIRPA = "Engulfment brake",
  C1QA = "Post-engulfment handling", C1QB = "Post-engulfment handling",
  C1QC = "Post-engulfment handling", APOE = "Post-engulfment handling",
  TREM2 = "Post-engulfment handling", CTSD = "Post-engulfment handling",
  LAMP1 = "Post-engulfment handling", NPC2 = "Post-engulfment handling",
  LIPA = "Post-engulfment handling", MAPK1 = "MEK–ERK transcript anchor",
  MAPK3 = "MEK–ERK transcript anchor", TGFB1 = "Repair and resolution",
  SMAD2 = "Repair and resolution", SMAD3 = "Repair and resolution",
  SMAD4 = "Repair and resolution"
)
candidate_plot <- unique(cross_data[
  article == "Article1" & cohort == "10_scRNA_discovery" &
    source_subset == "10_scRNA_old_discovery_DEG" &
    lineage %in% candidate_cell_order_plot & gene_norm %in% candidate_gene_order_plot,
  .(gene = gene_norm, cell_type = lineage, logFC = effect)
])
stopifnot(nrow(candidate_plot) == length(candidate_gene_order_plot) * length(candidate_cell_order_plot))
candidate_plot[, module := unname(candidate_module_plot[gene])]
candidate_module_colors <- c(
  "Hypoxia and entry" = set3[6], "Recognition and contact" = set3[5],
  "Engulfment brake" = set3[9], "Post-engulfment handling" = set3[1],
  "MEK–ERK transcript anchor" = set3[10], "Repair and resolution" = set3[8]
)
candidate_plot[, module_color := unname(candidate_module_colors[module])]
candidate_plot[, gene := factor(gene, levels = rev(candidate_gene_order_plot))]
candidate_plot[, cell_type := factor(cell_type, levels = candidate_cell_order_plot)]

p5b <- ggplot(candidate_plot, aes(x = logFC, y = gene, color = module)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  geom_segment(aes(x = 0, xend = logFC, yend = gene), linewidth = 0.55) +
  geom_point(size = 2.0) +
  facet_grid(. ~ cell_type) +
  scale_color_manual(values = candidate_module_colors) +
  scale_x_continuous(limits = c(-3.6, 1.5), breaks = c(-3, -2, -1, 0, 1)) +
  labs(
    title = "Discovery candidate-node contrasts",
    subtitle = "10-patient discovery cohort; severe versus mild logFC; color denotes functional grouping",
    x = "logFC", y = NULL, color = NULL
  ) +
  theme_pub +
  theme(
    axis.text.y = element_text(size = 6.3, face = "italic"),
    axis.text.x = element_text(size = 6.4),
    axis.ticks.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey96", linewidth = 0.2),
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 6.5, face = "bold", lineheight = 0.9),
    legend.position = "bottom", legend.key.width = unit(12, "mm"),
    legend.text = element_text(size = 5.9),
    plot.margin = margin(6, 7, 5, 7)
  )

# Panel C: a candidate relationship network rendered with points and arrows only.
axis_nodes <- data.table(
  x = c(0.75, 2.55, 4.35, 6.20, 8.05, 9.85), y = 0,
  label = c("HIF1A", "THBS1", "CD36", "MEK–ERK", "TGFB1", "SMAD2–4"),
  sublabel = c("", "", "", "MAPK1 and MAPK3\ntranscript anchor", "", ""),
  color = c(set3[6], set3[12], set3[5], set3[9], set3[1], set3[8])
)
arrow_data <- data.table(
  x = axis_nodes$x[-nrow(axis_nodes)] + 0.42,
  xend = axis_nodes$x[-1] - 0.42, y = 0, yend = 0,
  color = rep("grey55", nrow(axis_nodes) - 1)
)

p5c <- ggplot() +
  geom_segment(
    data = arrow_data,
    aes(x = x, xend = xend, y = y, yend = yend, color = color),
    linewidth = 0.55,
    linetype = "dashed",
    arrow = arrow(length = unit(1.5, "mm"), type = "open"),
    lineend = "round"
  ) +
  geom_segment(
    data = axis_nodes,
    aes(x = x - 0.38, xend = x + 0.38, y = -0.12, yend = -0.12, color = color),
    linewidth = 1.5, lineend = "round"
  ) +
  geom_text(
    data = axis_nodes,
    aes(x = x, y = 0.24, label = label, color = color),
    family = font_family, fontface = "bold", size = 2.6
  ) +
  geom_text(
    data = axis_nodes[nchar(sublabel) > 0],
    aes(x = x, y = -0.52, label = sublabel),
    family = font_family, size = 1.95, lineheight = 0.84, color = "grey30"
  ) +
  scale_color_identity() +
  coord_cartesian(xlim = c(0.10, 10.50), ylim = c(-1.15, 0.90), clip = "off") +
  labs(title = "HIF1A–THBS1–CD36–MEK–ERK–TGFB1–SMAD candidate network") +
  theme_void(base_family = font_family) +
  theme(
    plot.title = element_text(size = 9.4, face = "bold", hjust = 0),
    plot.margin = margin(6, 7, 4, 7)
  )

# Panel D: independent BALF proteome; negative mild-minus-severe means higher severe BALF abundance.
protein_targets <- c("HIF1A", "THBS1", "TGFB1", "AXL")
protein_plot <- unique(protein_data[
  comparison == "severity_Group1_vs_Group2" & gene_norm %in% protein_targets &
    !is.na(log2FC_group1_vs_group2),
  .(
    gene = gene_norm,
    effect = log2FC_group1_vs_group2,
    p_value,
    fdr = FDR
  )
])
stopifnot(nrow(protein_plot) == 4)
expected_protein <- c(HIF1A = -1.2401169, THBS1 = -0.3270706, TGFB1 = -0.1752031, AXL = 0.3975980)
stopifnot(all(abs(protein_plot$effect - expected_protein[protein_plot$gene]) < 1e-6))
protein_plot[, gene := factor(gene, levels = rev(protein_targets))]
protein_plot[, direction := fifelse(effect < 0, "Severe BALF higher", "Mild BALF higher")]
protein_plot[, stat_label := sprintf("P=%.4f; FDR=%.3f", p_value, fdr)]

p5d <- ggplot(protein_plot, aes(x = effect, y = gene, color = direction)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.45) +
  geom_segment(aes(x = 0, xend = effect, yend = gene), linewidth = 0.65) +
  geom_point(size = 2.6) +
  scale_color_manual(values = c("Severe BALF higher" = set3[4], "Mild BALF higher" = set3[5])) +
  scale_x_continuous(limits = c(-1.5, 1.18), breaks = c(-1.5, -1, -0.5, 0, 0.5, 1)) +
  labs(
    title = "Independent BALF proteome",
    subtitle = "110 BALF samples\nGroup1 mild 61; Group2 severe 49\nNominal P<0.05; FDR>0.05",
    x = "BALF protein log2 fold change\n(Group1 minus Group2)", y = NULL, color = NULL
  ) +
  theme_pub +
  guides(color = guide_legend(ncol = 1, byrow = TRUE)) +
  theme(
    legend.position = "bottom",
    legend.justification = "left",
    legend.text = element_text(size = 6.2),
    axis.title.x = element_text(size = 7.2, lineheight = 0.92),
    panel.grid.major.y = element_line(color = "grey96", linewidth = 0.2),
    plot.margin = margin(6, 9, 5, 7)
  )

# Panel E: concise point-based view of the parallel immune environment.
vdj_severity <- vdj_data[
  contrast == "重症_vs_轻症" & FDR_BH < 0.05,
  .(feature, log2_odds_ratio, FDR_BH)
]
stopifnot(nrow(vdj_severity) == 4)
stopifnot(all(c("TRBJ2-5", "TRAV3", "IGKV3D-20", "IGLV3-1") %in% vdj_severity$feature))
stopifnot(all(cellchat_data$fdr >= 0.05, na.rm = TRUE))

tcell_data <- data.table(
  feature = c("CSF1", "ITGA2", "MMP12"),
  logFC = c(0.717, -0.409, 2.005),
  fdr = c(0.864, 0.864, 0.864)
)
tcell_data[, feature := factor(feature, levels = c("MMP12", "ITGA2", "CSF1"))]
tcell_data[, `:=`(
  label_x = logFC + fifelse(logFC >= 0, 0.10, -0.10),
  label_hjust = fifelse(logFC >= 0, 0, 1)
)]
vdj_severity[, feature := as.character(feature)]
vdj_severity[, receptor := fifelse(feature %chin% c("TRBJ2-5", "TRAV3"), "TCR", "BCR")]
vdj_plot <- vdj_severity[match(c("TRBJ2-5", "TRAV3", "IGKV3D-20", "IGLV3-1"), feature)]
vdj_plot[, feature := factor(feature, levels = rev(c("TRBJ2-5", "TRAV3", "IGKV3D-20", "IGLV3-1")))]

mini_theme <- theme_classic(base_size = 6.6, base_family = font_family) +
  theme(
    plot.title = element_text(size = 7.2, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 5.8, color = "grey30", hjust = 0),
    axis.title = element_text(size = 6.1),
    axis.text = element_text(size = 5.8, color = "black"),
    axis.text.y = element_text(face = "italic"),
    axis.ticks.y = element_blank(),
    plot.margin = margin(2, 4, 2, 4)
  )

p5e_tcell <- ggplot(tcell_data, aes(x = logFC, y = feature)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.3) +
  geom_segment(aes(x = 0, xend = logFC, yend = feature), linewidth = 0.55, color = set3[5]) +
  geom_point(size = 2.0, color = set3[5]) +
  geom_text(aes(x = label_x, label = sprintf("%.3f", logFC), hjust = label_hjust),
            family = font_family, size = 1.65, color = "grey25") +
  scale_x_continuous(limits = c(-0.75, 2.45), breaks = c(0, 1, 2)) +
  labs(title = "T-cell pseudobulk", subtitle = "23-patient cohort; all three FDR=0.864", x = "logFC", y = NULL) +
  mini_theme

p5e_vdj <- ggplot(vdj_plot, aes(x = log2_odds_ratio, y = feature, color = receptor)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.3) +
  geom_segment(aes(x = 0, xend = log2_odds_ratio, yend = feature), linewidth = 0.55) +
  geom_point(size = 2.0) +
  scale_color_manual(values = c(TCR = set3[6], BCR = set3[1]), name = NULL) +
  scale_x_continuous(limits = c(-2.45, 0.75), breaks = c(-2, -1, 0)) +
  labs(title = "Severity-associated V(D)J", subtitle = "Four features at BH-FDR<0.05", x = "log2 odds ratio", y = NULL) +
  mini_theme +
  theme(legend.position = "bottom", legend.text = element_text(size = 5.6), legend.key.width = unit(7, "mm"))

context_data <- data.table(
  y = c(3, 2, 1),
  label = c(
    "Cytotoxic NKG7 down, FDR 0.039",
    "Neutrophil viral response down, FDR 0.0416",
    "CellChat: no pathway at BH-FDR<0.05"
  ),
  color = c(set3[5], set3[7], set3[9])
)
p5e_context <- ggplot(context_data, aes(x = 0, y = y)) +
  geom_point(aes(color = color), size = 1.9) +
  geom_text(aes(x = 0.08, label = label), hjust = 0, family = font_family, size = 1.92) +
  scale_color_identity() +
  coord_cartesian(xlim = c(-0.03, 1.24), ylim = c(0.55, 3.45), clip = "off") +
  labs(title = "Parallel immune context") +
  theme_void(base_family = font_family) +
  theme(
    plot.title = element_text(size = 7.2, face = "bold", hjust = 0),
    plot.margin = margin(2, 4, 2, 4)
  )

p5e <- wrap_elements(
  full = p5e_tcell / p5e_vdj / p5e_context + plot_layout(heights = c(1.0, 1.15, 0.78))
)

bottom_row <- p5d + p5e + plot_layout(widths = c(1.00, 1.00))

figure_5 <- p5a / p5b / p5c / bottom_row +
  plot_layout(heights = c(2.60, 3.60, 1.55, 3.10)) +
  plot_annotation(tag_levels = "A", theme = panel_tag_theme)

save_figure(figure_5, "Figure_5_cross_layer_candidate_model", width = 7.2, height = 10.3)

# Editable manuscript tables.
table_1 <- data.table(
  `数据资源` = c("10例发现单细胞队列", "23例独立扩大单细胞队列", "110例独立BALF蛋白组队列"),
  `轻症` = c(5L, 10L, 61L),
  `重症` = c(5L, 13L, 49L),
  `QC后分析规模` = c("115,594个细胞", "209,299个细胞", "110个独立生物学样本"),
  `主要作用` = c(
    "胞葬评分定位、候选节点和髓系状态排列",
    "独立扩大的免疫图谱、通路复现、组成、T细胞、CellChat和V(D)J",
    "HIF1A、THBS1、TGFB1及髓系稳态背景的跨层支持"
  ),
  `统计单位和配对关系` = c(
    "患儿；与其他队列独立",
    "患儿；不与10例合并为33例",
    "生物学样本；与单细胞患儿不做一对一配对"
  )
)
fwrite(table_1, file.path(table_dir, "Table_1_cohort_roles.csv"), bom = TRUE)

candidate_genes <- c(
  "HIF1A", "THBS1", "CD36", "MERTK", "AXL", "CD47", "SIRPA",
  "C1QA", "C1QB", "C1QC", "APOE", "TREM2", "CTSD", "LAMP1", "NPC2", "LIPA",
  "MAPK1", "MAPK3", "TGFB1", "SMAD2", "SMAD3", "SMAD4"
)
gene_modules <- c(
  HIF1A = "缺氧与代谢入口",
  THBS1 = "识别、接触和吞入", CD36 = "识别、接触和吞入",
  MERTK = "识别、接触和吞入", AXL = "识别、接触和吞入",
  CD47 = "吞入制动", SIRPA = "吞入制动",
  C1QA = "吞入后处理", C1QB = "吞入后处理", C1QC = "吞入后处理",
  APOE = "吞入后处理", TREM2 = "吞入后处理", CTSD = "吞入后处理",
  LAMP1 = "吞入后处理", NPC2 = "吞入后处理", LIPA = "吞入后处理",
  MAPK1 = "MEK–ERK候选中继转录节点", MAPK3 = "MEK–ERK候选中继转录节点",
  TGFB1 = "修复和炎症收束", SMAD2 = "修复和炎症收束",
  SMAD3 = "修复和炎症收束", SMAD4 = "修复和炎症收束"
)
target_cell_types <- c("Mono-derived Macro", "Mono-derived alveolar Macro", "Monocytes")
candidate_long <- unique(cross_data[
  article == "Article1" & cohort == "10_scRNA_discovery" &
    source_subset == "10_scRNA_old_discovery_DEG" &
    lineage %in% target_cell_types & gene_norm %in% candidate_genes,
  .(gene = gene_norm, cell_type = lineage, logFC = effect)
])
stopifnot(nrow(candidate_long) == length(candidate_genes) * length(target_cell_types))
candidate_wide <- dcast(candidate_long, gene ~ cell_type, value.var = "logFC")
setnames(
  candidate_wide,
  target_cell_types,
  c("单核来源巨噬细胞logFC", "单核来源肺泡样巨噬细胞logFC", "单核细胞logFC")
)
candidate_wide[, `功能模块` := unname(gene_modules[gene])]
candidate_wide[, `节点` := gene]
candidate_wide[, gene := NULL]
setcolorder(candidate_wide, c(
  "功能模块", "节点", "单核来源巨噬细胞logFC",
  "单核来源肺泡样巨噬细胞logFC", "单核细胞logFC"
))
candidate_wide[, (names(candidate_wide)[3:5]) := lapply(.SD, round, 3), .SDcols = 3:5]
candidate_wide[, `节点` := factor(`节点`, levels = candidate_genes)]
setorder(candidate_wide, `节点`)
candidate_wide[, `节点` := as.character(`节点`)]
fwrite(candidate_wide, file.path(table_dir, "Table_2_candidate_nodes.csv"), bom = TRUE)

table_3 <- copy(protein_plot)
table_3[, gene := as.character(gene)]
table_3[, `候选节点` := gene]
table_3[, `Group1轻症-Group2重症 log2FC` := round(effect, 3)]
table_3[, `名义P值` := signif(p_value, 3)]
table_3[, `FDR` := round(fdr, 3)]
table_3[, `方向` := fifelse(effect < 0, "重症BALF整体蛋白相对更高", "重症BALF整体蛋白相对更低")]
table_3 <- table_3[, .(
  `候选节点`, `Group1轻症-Group2重症 log2FC`, `名义P值`, FDR, `方向`
)]
table_3[, sort_order := match(`候选节点`, protein_targets)]
setorder(table_3, sort_order)
table_3[, sort_order := NULL]
fwrite(table_3, file.path(table_dir, "Table_3_proteomic_direction.csv"), bom = TRUE)

traceability <- data.table(
  artifact = c(
    "Figure 5A", "Figure 5B-C", "Figure 5D", "Figure 5E T-cell candidates",
    "Figure 5E CellChat", "Figure 5E V(D)J", "Table 1", "Table 2", "Table 3"
  ),
  source = c(
    axis_file, cross_file, protein_file,
    file.path(source_dir, "两篇文章定稿前深度审计报告_2026-08-20.md"),
    cellchat_file, vdj_file,
    file.path(source_dir, "两篇文章_方法与统计口径补全_2026-08-19.md"),
    cross_file, protein_file
  ),
  scope = c(
    "10-patient discovery score means",
    "Discovery candidate-node organization",
    "110-sample independent BALF protein effects",
    "Locked 23-patient T-cell results",
    "Patient-level candidate interaction pathways",
    "Patient-level severity-associated receptor features",
    "Locked cohort identities and analysis roles",
    "Discovery candidate-node logFC values",
    "Independent BALF protein effects and direction"
  )
)
fwrite(traceability, file.path(table_dir, "Figure_table_source_traceability.csv"), bom = TRUE)

cat("Built Figure 5 and manuscript tables in:\n", output_dir, "\n")
