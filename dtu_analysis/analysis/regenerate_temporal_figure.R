#!/usr/bin/env Rscript

# Regenerate Figure 2 from the committed contrast-count table without rerunning
# the archive-level analysis. The complete workflow uses the same plotting code
# in run_analysis.R.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- sub("^--file=", "", script_argument[[1]])
paper_dir <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/")
figure_dir <- file.path(paper_dir, "figures")
table_dir <- file.path(paper_dir, "tables")

tissue_colors <- c(Forebrain = "#0072B2", Hindbrain = "#009E73", Midbrain = "#D55E00")
tissue_shapes <- c(Forebrain = 16, Hindbrain = 17, Midbrain = 15)
tissue_linetypes <- c(Forebrain = "solid", Hindbrain = "22", Midbrain = "42")
stage_order <- c("10.5", "11.5", "12.5", "13.5", "14.5", "15.5", "16.5", "0")
stage_labels <- c(
  "10.5" = "E10.5", "11.5" = "E11.5", "12.5" = "E12.5", "13.5" = "E13.5",
  "14.5" = "E14.5", "15.5" = "E15.5", "16.5" = "E16.5", "0" = "P0"
)

theme_paper <- function(base_size = 10) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = rel(1.05)),
      plot.subtitle = element_text(color = "#4D4D4D"),
      axis.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      plot.margin = margin(5.5, 8, 5.5, 5.5)
    )
}

contrast_counts <- read.csv(file.path(table_dir, "temporal_contrast_counts.csv"))
contrast_counts <- contrast_counts |>
  mutate(
    early = as.character(early),
    late = as.character(late),
    tissue = factor(tissue, levels = names(tissue_colors))
  )

p_heat <- contrast_counts |>
  mutate(
    early = factor(early, levels = stage_order, labels = unname(stage_labels[stage_order])),
    late = factor(late, levels = stage_order, labels = unname(stage_labels[stage_order]))
  ) |>
  ggplot(aes(late, early, fill = genes)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = genes), size = 2.55, color = "#1A1A1A") +
  facet_wrap(~ tissue, nrow = 1) +
  scale_fill_gradient(low = "#F7FBFF", high = "#08519C", name = "DTU genes") +
  coord_equal() +
  labs(
    title = "A  DTU counts across developmental stage contrasts",
    subtitle = "Unique significant genes in each of the 28 pairwise stage contrasts",
    x = "Later stage", y = "Earlier stage"
  ) +
  theme_paper(11) +
  theme(axis.text.x = element_text(angle = 55, hjust = 1), legend.position = "bottom")

adjacent <- contrast_counts |>
  filter(late_i - early_i == 1) |>
  mutate(
    interval = factor(
      paste(stage_labels[early], stage_labels[late], sep = "\u2013"),
      levels = paste(
        stage_labels[stage_order[-length(stage_order)]],
        stage_labels[stage_order[-1]],
        sep = "\u2013"
      )
    )
  )

p_adjacent <- adjacent |>
  ggplot(aes(interval, genes, color = tissue, shape = tissue, linetype = tissue, group = tissue)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_color_manual(values = tissue_colors, name = NULL) +
  scale_shape_manual(values = tissue_shapes, name = NULL) +
  scale_linetype_manual(values = tissue_linetypes, name = NULL) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "B  Adjacent-stage contrasts expose an E15.5-centred midbrain peak",
    x = "Developmental interval", y = "DTU genes"
  ) +
  theme_paper(11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")

figure2 <- p_heat / p_adjacent + plot_layout(heights = c(1.55, 0.8))
ggsave(
  file.path(figure_dir, "figure2_temporal_dtu.pdf"), figure2,
  width = 9.2, height = 7.0, units = "in", device = cairo_pdf
)
ggsave(
  file.path(figure_dir, "figure2_temporal_dtu.png"), figure2,
  width = 9.2, height = 7.0, units = "in", dpi = 300, bg = "white"
)
