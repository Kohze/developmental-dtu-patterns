#!/usr/bin/env Rscript

# Regenerate the trajectory-only main figure from the committed candidate table.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else "."
paper_dir <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/")
table_path <- file.path(
  paper_dir,
  "tables",
  "transient_regional_top_candidate_trajectories.csv"
)
figure_dir <- file.path(paper_dir, "figures")

stage_levels <- c("10.5", "11.5", "12.5", "13.5", "14.5", "15.5", "16.5", "0")
stage_labels <- c(
  `10.5` = "E10.5", `11.5` = "E11.5", `12.5` = "E12.5",
  `13.5` = "E13.5", `14.5` = "E14.5", `15.5` = "E15.5",
  `16.5` = "E16.5", `0` = "P0"
)

candidate_means <- read.csv(table_path, stringsAsFactors = FALSE) |>
  mutate(stage = factor(stage, levels = stage_levels))

facet_labels <- c(
  Ntrk2 = "Ntrk2\nhigher NM_008745\nlower NM_001025074",
  Scg3 = "Scg3\nhigher NM_009130\nlower NM_001164790",
  Tecr = "Tecr\nhigher NM_134118\nlower NM_027179",
  Armc8 = "Armc8\nhigher NM_028768\nlower NM_001166138",
  Bin1 = "Bin1\nhigher NM_001083334\nlower NM_009668",
  Gpm6a = "Gpm6a\nhigher NM_001253754\nlower NM_153581"
)

figure <- ggplot(
  candidate_means,
  aes(
    stage,
    isoform_fraction,
    colour = direction,
    linetype = region_group,
    group = interaction(isoform_id, region_group)
  )
) +
  annotate(
    "rect", xmin = 5.5, xmax = 6.5, ymin = -Inf, ymax = Inf,
    fill = "#F0E442", alpha = 0.12
  ) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.35) +
  facet_wrap(
    ~ gene_name, scales = "free_y", ncol = 3,
    labeller = as_labeller(facet_labels)
  ) +
  scale_colour_manual(
    values = c(higher = "#D55E00", lower = "#0072B2"),
    labels = c(higher = "Higher at episode", lower = "Lower at episode")
  ) +
  scale_x_discrete(labels = stage_labels) +
  scale_linetype_manual(
    values = c("Midbrain" = "solid", "Forebrain/hindbrain mean" = "22")
  ) +
  labs(
    title = "Leading replicate-consistent reciprocal exchanges",
    subtitle = "Higher/lower denote the midbrain fraction at E15.5; RefSeq accessions are shown in each facet",
    x = "Developmental stage",
    y = "Isoform fraction",
    colour = NULL,
    linetype = NULL
  ) +
  theme_minimal(base_size = 9.5) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(figure_dir, "figure6b_top_candidate_trajectories.pdf"),
  figure,
  width = 7.2,
  height = 5.4,
  units = "in",
  device = cairo_pdf
)
ggsave(
  file.path(figure_dir, "figure6b_top_candidate_trajectories.png"),
  figure,
  width = 7.2,
  height = 5.4,
  units = "in",
  dpi = 300,
  bg = "white"
)
