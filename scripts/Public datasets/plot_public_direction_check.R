#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

base_dir <- Sys.getenv("MPP_PUBLIC_ROOT", "inputs/public")
result_dir <- file.path(base_dir, "results")
figure_dir <- file.path(base_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

gse274 <- fread(file.path(result_dir, "GSE274269_candidate_nodes.csv"))
gse179 <- fread(file.path(result_dir, "GSE179051_candidate_nodes.csv"))

node_order <- gse274$gene
gse274[, gene := factor(gene, levels = rev(node_order))]
gse179[, gene := factor(gene, levels = rev(node_order))]

set3 <- c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5")

p_detection <- ggplot(gse274, aes(x = gene, y = detected_rate * 100, fill = "BALF total RNA")) +
  geom_col(width = 0.68, colour = "black", linewidth = 0.18) +
  coord_flip(clip = "off") +
  scale_fill_manual(values = set3[1], name = NULL) +
  scale_y_continuous(limits = c(0, 105), breaks = c(0, 25, 50, 75, 100),
                     labels = function(x) paste0(x, "%"), expand = c(0, 0)) +
  labs(x = "Candidate node", y = "Detection rate") +
  theme_classic(base_family = "Times New Roman", base_size = 10.5) +
  theme(legend.position = c(0.72, 0.08), legend.background = element_blank(),
        axis.text.y = element_text(colour = "black", size = 8.6),
        axis.text.x = element_text(colour = "black"), axis.title = element_text(colour = "black"),
        plot.margin = margin(5.5, 10, 5.5, 5.5))

long_fc <- rbindlist(list(
  gse274[, .(gene, logFC = logFC_P1_2_minus_P1_1,
             dataset = "GSE274269 BALF; P1-2 minus P1-1")],
  gse179[, .(gene, logFC = logFC_M_pneumoniae_minus_PBS,
             dataset = "GSE179051 HeLa; M. pneumoniae minus PBS")]
), fill = TRUE)
long_fc[, dataset := factor(dataset, levels = c(
  "GSE274269 BALF; P1-2 minus P1-1",
  "GSE179051 HeLa; M. pneumoniae minus PBS"
))]

p_direction <- ggplot(long_fc, aes(x = logFC, y = gene, colour = dataset)) +
  geom_vline(xintercept = 0, colour = "black", linewidth = 0.3) +
  geom_point(size = 2.35, na.rm = TRUE) +
  scale_colour_manual(values = set3[c(3, 5)], name = NULL) +
  labs(x = "logFC", y = NULL) +
  theme_classic(base_family = "Times New Roman", base_size = 10.5) +
  theme(legend.position = "bottom", legend.background = element_blank(),
        axis.text.y = element_text(colour = "black", size = 8.6),
        axis.text.x = element_text(colour = "black"), axis.title = element_text(colour = "black"),
        plot.margin = margin(5.5, 10, 5.5, 5.5))

figure <- p_detection + p_direction + plot_layout(widths = c(1, 1.05))
figure <- figure + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(family = "Times New Roman", face = "bold", size = 12,
                                colour = "black"), plot.tag.position = c(0.005, 0.995))

stem <- file.path(figure_dir, "Supplementary_Figure_public_direction_check")
ggsave(paste0(stem, ".svg"), figure, width = 11.2, height = 7.1, units = "in", device = "svg")
ggsave(paste0(stem, ".png"), figure, width = 11.2, height = 7.1, units = "in", dpi = 300, bg = "white")
ggsave(paste0(stem, ".tiff"), figure, width = 11.2, height = 7.1, units = "in", dpi = 600, compression = "lzw", bg = "white")
wrapper_html <- c(
  "<!doctype html>",
  "<html>",
  "<head>",
  "<meta charset='utf-8'>",
  "<style>",
  "  @page { size: 11.2in 7.1in; margin: 0; }",
  "  html, body { width: 11.2in; height: 7.1in; margin: 0; padding: 0; background: #fff; }",
  "  img { display: block; width: 11.2in; height: 7.1in; object-fit: contain; }",
  "</style>",
  "</head>",
  "<body>",
  "  <img src='Supplementary_Figure_public_direction_check.svg' alt='Public database direction check figure'>",
  "</body>",
  "</html>"
)
writeLines(wrapper_html, file.path(figure_dir, "public_direction_check_pdf_wrapper.html"))
render_status <- system2(
  "python3",
  shQuote(c(file.path(base_dir, "render_public_direction_check_pdf.py"),
            paste0(stem, ".svg"),
            paste0(stem, ".pdf"),
            file.path(figure_dir, "public_direction_check_pdf_wrapper.html"))),
  stdout = TRUE, stderr = TRUE
)
cat(render_status, sep = "\n")
writeLines(capture.output(sessionInfo()), file.path(figure_dir, "Supplementary_Figure_public_direction_check_sessionInfo.txt"))
cat("Wrote public direction-check figure files to", figure_dir, "\n")
