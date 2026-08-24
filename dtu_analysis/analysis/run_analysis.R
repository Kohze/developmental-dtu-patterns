#!/usr/bin/env Rscript

# Reproducible analysis for the DTU manuscript.
# Usage:
#   Rscript analysis/run_analysis.R [directory containing the two .RData inputs]

options(stringsAsFactors = FALSE, width = 120)
set.seed(20260729)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(patchwork)
  library(clusterProfiler)
  library(org.Mm.eg.db)
})

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_argument)) {
  stop("Cannot determine the script path. Run this file with Rscript.")
}
script_file <- sub("^--file=", "", script_argument[[1]])
paper_dir <- normalizePath(
  file.path(dirname(script_file), ".."),
  winslash = "/",
  mustWork = TRUE
)

figure_dir <- file.path(paper_dir, "figures")
table_dir <- file.path(paper_dir, "tables")
data_dir <- file.path(paper_dir, "data")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

required_inputs <- c("foreBrain.RData", "isoform_final.RData")
cli_arguments <- commandArgs(trailingOnly = TRUE)
if (length(cli_arguments) > 1L) {
  stop("Expected at most one argument: the directory containing the .RData inputs.")
}
input_candidates <- if (length(cli_arguments)) {
  cli_arguments[[1]]
} else {
  c(paper_dir, file.path(paper_dir, "..", ".."), getwd())
}
input_candidates <- unique(normalizePath(
  input_candidates,
  winslash = "/",
  mustWork = FALSE
))
complete_input <- vapply(
  input_candidates,
  function(path) all(file.exists(file.path(path, required_inputs))),
  logical(1)
)
if (!any(complete_input)) {
  stop(
    "Could not find both required inputs (",
    paste(required_inputs, collapse = ", "),
    "). Pass their directory as the sole command-line argument."
  )
}
input_dir <- normalizePath(
  input_candidates[which(complete_input)[[1]]],
  winslash = "/",
  mustWork = TRUE
)
message("Loading archived inputs from: ", input_dir)

load(file.path(input_dir, "foreBrain.RData"))
load(file.path(input_dir, "isoform_final.RData"))

objects <- list(
  Forebrain = combinedForebrain,
  Hindbrain = combinedHindBrain,
  Midbrain = combinedMitbrain
)
background_objects <- list(
  Forebrain = switchListForebrain,
  Hindbrain = switchListHindbrain,
  Midbrain = switchListMidbrain
)

tissue_colors <- c(
  Forebrain = "#0072B2",
  Hindbrain = "#009E73",
  Midbrain = "#D55E00"
)
tissue_shapes <- c(
  Forebrain = 16,
  Hindbrain = 17,
  Midbrain = 15
)
tissue_linetypes <- c(
  Forebrain = "solid",
  Hindbrain = "22",
  Midbrain = "42"
)
stage_order <- c("10.5", "11.5", "12.5", "13.5", "14.5", "15.5", "16.5", "0")
stage_labels <- c(
  "10.5" = "E10.5", "11.5" = "E11.5", "12.5" = "E12.5",
  "13.5" = "E13.5", "14.5" = "E14.5", "15.5" = "E15.5",
  "16.5" = "E16.5", "0" = "P0"
)
event_labels <- c(
  ES = "Exon skipping", MEE = "Mutually exclusive exon",
  MES = "Multiple exon skipping", IR = "Intron retention",
  A5 = "Alternative 5' splice site", A3 = "Alternative 3' splice site",
  ATSS = "Alternative transcription start", ATTS = "Alternative transcription termination"
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

save_figure <- function(plot, stem, width, height) {
  ggsave(file.path(figure_dir, paste0(stem, ".pdf")), plot,
         width = width, height = height, units = "in", device = cairo_pdf)
  ggsave(file.path(figure_dir, paste0(stem, ".png")), plot,
         width = width, height = height, units = "in", dpi = 300, bg = "white")
}

is_significant <- function(d) {
  !is.na(d$isoform_switch_q_value) &
    d$isoform_switch_q_value < 0.05 &
    abs(d$dIF) >= 0.10
}

sig_rows <- lapply(objects, function(x) {
  d <- as.data.frame(x$isoformFeatures)
  d[is_significant(d), , drop = FALSE]
})
gene_sets <- lapply(sig_rows, function(d) unique(d$gene_id))
isoform_sets <- lapply(sig_rows, function(d) unique(d$isoform_id))
gene_union <- Reduce(union, gene_sets)
gene_core <- Reduce(intersect, gene_sets)
midbrain_only <- setdiff(gene_sets$Midbrain, union(gene_sets$Forebrain, gene_sets$Hindbrain))
analysis_universe <- Reduce(
  union,
  lapply(background_objects, function(x) unique(x$isoformFeatures$gene_id))
)

# Dataset and DTU summary -----------------------------------------------------
dataset_summary <- bind_rows(lapply(names(objects), function(tissue) {
  x <- objects[[tissue]]
  data.frame(
    Tissue = tissue,
    Libraries = nrow(x$designMatrix),
    Stages = length(unique(x$designMatrix$condition)),
    DTU_genes = length(gene_sets[[tissue]]),
    DTU_isoforms = length(isoform_sets[[tissue]]),
    Significant_isoform_contrasts = nrow(sig_rows[[tissue]])
  )
}))
write.csv(dataset_summary, file.path(table_dir, "dataset_summary.csv"), row.names = FALSE)

membership <- data.frame(
  gene_id = gene_union,
  Forebrain = gene_union %in% gene_sets$Forebrain,
  Hindbrain = gene_union %in% gene_sets$Hindbrain,
  Midbrain = gene_union %in% gene_sets$Midbrain
) |>
  mutate(
    membership = case_when(
      Forebrain & Hindbrain & Midbrain ~ "Forebrain + Hindbrain + Midbrain",
      Forebrain & Hindbrain ~ "Forebrain + Hindbrain",
      Forebrain & Midbrain ~ "Forebrain + Midbrain",
      Hindbrain & Midbrain ~ "Hindbrain + Midbrain",
      Forebrain ~ "Forebrain only",
      Hindbrain ~ "Hindbrain only",
      Midbrain ~ "Midbrain only"
    )
  )
write.csv(membership, file.path(data_dir, "dtu_gene_membership.csv"), row.names = FALSE)

membership_counts <- membership |>
  dplyr::count(membership, name = "genes") |>
  arrange(desc(genes)) |>
  mutate(membership = factor(membership, levels = rev(membership)))
write.csv(membership_counts, file.path(table_dir, "dtu_membership_counts.csv"), row.names = FALSE)

p_region <- dataset_summary |>
  mutate(Tissue = factor(Tissue, levels = rev(names(tissue_colors)))) |>
  ggplot(aes(DTU_genes, Tissue, fill = Tissue)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  geom_text(aes(label = format(DTU_genes, big.mark = ",")), hjust = -0.12, size = 3.4) +
  scale_fill_manual(values = tissue_colors) +
  scale_x_continuous(limits = c(0, 2700), expand = expansion(mult = c(0, 0))) +
  labs(
    title = "A  DTU genes by region",
    subtitle = "Unique genes significant in at least one temporal contrast",
    x = "DTU genes", y = NULL
  ) +
  theme_paper()

p_membership <- ggplot(membership_counts, aes(genes, membership)) +
  geom_col(fill = "#5B6470", width = 0.65) +
  geom_text(aes(label = format(genes, big.mark = ",")), hjust = -0.12, size = 3.4) +
  scale_x_continuous(limits = c(0, 2200), expand = expansion(mult = c(0, 0))) +
  labs(
    title = "B  A shared core and a midbrain detection-only component",
    subtitle = paste0(format(length(gene_union), big.mark = ","),
                      " DTU genes in the three-region union"),
    x = "Genes", y = NULL
  ) +
  theme_paper()

figure1 <- p_region / p_membership + plot_layout(heights = c(0.8, 1.35))
save_figure(figure1, "figure1_dtu_landscape", 7.2, 6.2)

# Temporal contrast landscape ------------------------------------------------
clean_condition <- function(x) sub("^X", "", x)

contrast_counts <- bind_rows(lapply(names(sig_rows), function(tissue) {
  d <- sig_rows[[tissue]] |>
    transmute(
      tissue = tissue,
      gene_id,
      c1 = clean_condition(condition_1),
      c2 = clean_condition(condition_2)
    ) |>
    mutate(
      i1 = match(c1, stage_order),
      i2 = match(c2, stage_order),
      early_i = pmin(i1, i2),
      late_i = pmax(i1, i2),
      early = stage_order[early_i],
      late = stage_order[late_i]
    ) |>
    distinct(tissue, gene_id, early_i, late_i, early, late) |>
    dplyr::count(tissue, early_i, late_i, early, late, name = "genes")
  d
}))
write.csv(contrast_counts, file.path(table_dir, "temporal_contrast_counts.csv"), row.names = FALSE)

p_heat <- contrast_counts |>
  mutate(
    early = factor(early, levels = stage_order, labels = unname(stage_labels[stage_order])),
    late = factor(late, levels = stage_order, labels = unname(stage_labels[stage_order])),
    tissue = factor(tissue, levels = names(tissue_colors))
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
      levels = paste(stage_labels[stage_order[-length(stage_order)]],
                     stage_labels[stage_order[-1]], sep = "\u2013")
    )
  )
write.csv(adjacent, file.path(table_dir, "adjacent_stage_counts.csv"), row.names = FALSE)

p_adjacent <- adjacent |>
  ggplot(aes(
    interval, genes,
    color = tissue, shape = tissue, linetype = tissue, group = tissue
  )) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_color_manual(values = tissue_colors, name = NULL) +
  scale_shape_manual(values = tissue_shapes, name = NULL) +
  scale_linetype_manual(values = tissue_linetypes, name = NULL) +
  scale_y_continuous(
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    title = "B  Adjacent-stage contrasts expose an E15.5-centred midbrain peak",
    x = "Developmental interval", y = "DTU genes"
  ) +
  theme_paper(11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )

figure2 <- p_heat / p_adjacent + plot_layout(heights = c(1.55, 0.8))
save_figure(figure2, "figure2_temporal_dtu", 9.2, 7.0)

# Transcript-structure event composition ------------------------------------
event_cols <- names(event_labels)
event_prevalence <- bind_rows(lapply(names(objects), function(tissue) {
  a <- as.data.frame(objects[[tissue]]$AlternativeSplicingAnalysis)
  a <- a[a$isoform_id %in% isoform_sets[[tissue]], c("isoform_id", event_cols)]
  bind_rows(lapply(event_cols, function(event) {
    value <- a[[event]]
    present <- if (is.logical(value)) value else !is.na(value) & value != 0
    data.frame(
      tissue = tissue,
      event = event,
      isoforms = sum(present, na.rm = TRUE),
      total_isoforms = nrow(a),
      prevalence = mean(present, na.rm = TRUE)
    )
  }))
}))
write.csv(event_prevalence, file.path(table_dir, "splicing_event_prevalence.csv"), row.names = FALSE)

figure3 <- event_prevalence |>
  mutate(
    event = factor(event_labels[event], levels = rev(unname(event_labels))),
    tissue = factor(tissue, levels = names(tissue_colors))
  ) |>
  ggplot(aes(prevalence, event, color = tissue, shape = tissue)) +
  geom_line(aes(group = event), color = "#D9D9D9", linewidth = 0.6) +
  geom_point(size = 2.7, position = position_dodge(width = 0.35)) +
  scale_color_manual(values = tissue_colors, name = NULL) +
  scale_shape_manual(values = tissue_shapes, name = NULL) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, max(event_prevalence$prevalence) * 1.08)) +
  labs(
    title = "DTU isoforms carry multiple transcript-structure annotations",
    subtitle = "Events are non-exclusive; denominator is the significant isoform set in each region",
    x = "Fraction of DTU isoforms carrying event", y = NULL
  ) +
  theme_paper(9) +
  theme(legend.position = "top")
save_figure(figure3, "figure3_splicing_architecture", 7.2, 4.8)

# Functional differentiation within the DTU universe ------------------------
go_sets <- list(
  `Shared core` = gene_core,
  `Midbrain only` = midbrain_only
)

go_results <- bind_rows(lapply(names(go_sets), function(group) {
  genes <- go_sets[[group]]
  ego <- enrichGO(
    gene = genes,
    universe = analysis_universe,
    OrgDb = org.Mm.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.20,
    minGSSize = 10,
    maxGSSize = 500,
    readable = FALSE
  )
  result <- as.data.frame(ego)
  if (!nrow(result)) return(data.frame())
  result$group <- group
  result
}))

if (nrow(go_results)) {
  parse_ratio <- function(x) {
    bits <- strsplit(x, "/", fixed = TRUE)
    vapply(bits, function(z) as.numeric(z[1]) / as.numeric(z[2]), numeric(1))
  }
  go_results <- go_results |>
    mutate(
      GeneFraction = parse_ratio(GeneRatio),
      BackgroundFraction = parse_ratio(BgRatio),
      enrichment = GeneFraction / BackgroundFraction
    )
  write.csv(go_results, file.path(data_dir, "go_bp_enrichment_all.csv"), row.names = FALSE)

  top_go <- go_results |>
    group_by(group) |>
    arrange(p.adjust, desc(enrichment), .by_group = TRUE) |>
    slice_head(n = 8) |>
    ungroup() |>
    mutate(
      term = str_wrap(Description, width = 42),
      term = factor(term, levels = rev(unique(term))),
      group = factor(group, levels = c("Shared core", "Midbrain only"))
    )
  write.csv(top_go, file.path(table_dir, "go_bp_enrichment_top.csv"), row.names = FALSE)

  figure4 <- ggplot(top_go, aes(enrichment, term, size = Count, color = -log10(p.adjust))) +
    geom_point(alpha = 0.9) +
    facet_wrap(~ group, scales = "free_y", nrow = 1) +
    scale_color_viridis_c(option = "C", end = 0.88, name = expression(-log[10]("FDR"))) +
    scale_size_continuous(range = c(2.2, 7), name = "Genes") +
    labs(
      title = "DTU genes are enriched for neuronal morphogenesis and plasticity",
      subtitle = sprintf(
        "GO biological-process over-representation against %s background genes",
        format(length(analysis_universe), big.mark = ",")
      ),
      x = "Fold enrichment", y = NULL
    ) +
    theme_paper(10.5) +
    theme(legend.position = "bottom")
  save_figure(figure4, "figure4_go_enrichment", 9.2, 6.0)
}

# Representative midbrain isoform-usage trajectories ------------------------
isoform_gene_map <- objects$Midbrain$isoformFeatures |>
  as.data.frame() |>
  distinct(gene_id, gene_name, isoform_id)

make_trajectory <- function(gene) {
  ids <- isoform_gene_map |>
    filter(gene_id == gene) |>
    pull(isoform_id) |>
    unique()
  d <- as.data.frame(objects$Midbrain$isoformRepIF)
  d <- d[d$isoform_id %in% ids, , drop = FALSE]
  if (!nrow(d)) return(data.frame())
  long <- d |>
    pivot_longer(-isoform_id, names_to = "sample", values_to = "IF") |>
    separate(sample, into = c("stage", "replicate"), sep = "_", remove = FALSE) |>
    group_by(isoform_id) |>
    mutate(mean_overall = mean(IF, na.rm = TRUE)) |>
    ungroup()
  keep <- long |>
    distinct(isoform_id, mean_overall) |>
    arrange(desc(mean_overall)) |>
    slice_head(n = 3) |>
    pull(isoform_id)
  long |>
    filter(isoform_id %in% keep) |>
    mutate(
      gene = gene,
      stage = factor(stage, levels = stage_order,
                     labels = unname(stage_labels[stage_order]))
    ) |>
    group_by(gene, isoform_id, stage) |>
    summarise(
      mean_IF = mean(IF, na.rm = TRUE),
      min_IF = min(IF, na.rm = TRUE),
      max_IF = max(IF, na.rm = TRUE),
      .groups = "drop"
    )
}

trajectory_genes <- c("Ppp2r3a", "Rtn2")
trajectory <- bind_rows(lapply(trajectory_genes, make_trajectory))
write.csv(trajectory, file.path(table_dir, "representative_isoform_trajectories.csv"), row.names = FALSE)

if (nrow(trajectory)) {
  isoform_colors <- c(
    "NM_001161362" = "#0072B2",
    "NM_172144" = "#D55E00",
    "NM_001025364" = "#009E73",
    "NM_013648" = "#CC79A7"
  )
  isoform_shapes <- c(
    "NM_001161362" = 16,
    "NM_172144" = 17,
    "NM_001025364" = 15,
    "NM_013648" = 18
  )
  isoform_linetypes <- c(
    "NM_001161362" = "solid",
    "NM_172144" = "22",
    "NM_001025364" = "42",
    "NM_013648" = "44"
  )
  figure5 <- ggplot(
    trajectory,
    aes(
      stage, mean_IF,
      color = isoform_id, shape = isoform_id,
      linetype = isoform_id, group = isoform_id
    )
  ) +
    geom_ribbon(aes(ymin = min_IF, ymax = max_IF, fill = isoform_id),
                alpha = 0.12, color = NA, show.legend = FALSE) +
    geom_line(linewidth = 0.85) +
    geom_point(size = 1.8) +
    facet_wrap(~ gene, nrow = 1) +
    scale_color_manual(
      values = isoform_colors,
      limits = names(isoform_colors),
      name = "Archived RefSeq accession"
    ) +
    scale_fill_manual(values = isoform_colors, limits = names(isoform_colors)) +
    scale_shape_manual(
      values = isoform_shapes,
      limits = names(isoform_shapes),
      name = "Archived RefSeq accession"
    ) +
    scale_linetype_manual(
      values = isoform_linetypes,
      limits = names(isoform_linetypes),
      name = "Archived RefSeq accession"
    ) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       limits = c(0, 1)) +
    labs(
      title = "Original heuristic candidates show transient E15.5 redistribution",
      subtitle = "Historical comparison; lines are replicate means and ribbons span both replicates",
      x = "Developmental stage", y = "Isoform fraction"
    ) +
    theme_paper(10.5) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(face = "bold.italic"),
      legend.position = "bottom"
    )
  save_figure(figure5, "figure5_candidate_trajectories", 8.3, 4.9)
}

# Sample-level audit ----------------------------------------------------------
expression_matrix <- function(x) {
  d <- as.data.frame(x$isoformRepExpression)
  rownames(d) <- d$isoform_id
  as.matrix(d[, setdiff(names(d), "isoform_id"), drop = FALSE])
}
expression_mats <- lapply(background_objects, expression_matrix)
common_isoforms <- Reduce(intersect, lapply(expression_mats, rownames))

correlation_audit <- bind_rows(lapply(names(expression_mats), function(tissue) {
  m <- expression_mats[[tissue]][common_isoforms, , drop = FALSE]
  bind_rows(lapply(stage_order, function(stage) {
    data.frame(
      comparison = "Within-region biological replicates",
      region_pair = tissue,
      stage = stage_labels[stage],
      replicate = "1 versus 2",
      pearson_log1p = cor(log1p(m[, paste0(stage, "_1")]),
                         log1p(m[, paste0(stage, "_2")]),
                         use = "pairwise.complete.obs")
    )
  }))
}), bind_rows(lapply(list(c("Forebrain", "Hindbrain"),
                         c("Forebrain", "Midbrain"),
                         c("Hindbrain", "Midbrain")), function(pair) {
  a <- expression_mats[[pair[1]]][common_isoforms, , drop = FALSE]
  b <- expression_mats[[pair[2]]][common_isoforms, , drop = FALSE]
  bind_rows(lapply(stage_order, function(stage) {
    bind_rows(lapply(1:2, function(rep_a) {
      bind_rows(lapply(1:2, function(rep_b) {
        sample_a <- paste0(stage, "_", rep_a)
        sample_b <- paste0(stage, "_", rep_b)
        data.frame(
          comparison = "Same-stage cross-region replicate combinations",
          region_pair = paste(pair, collapse = " vs "),
          stage = stage_labels[stage],
          replicate = paste(rep_a, "versus", rep_b),
          pearson_log1p = cor(log1p(a[, sample_a]), log1p(b[, sample_b]),
                             use = "pairwise.complete.obs")
        )
      }))
    }))
  }))
})))
write.csv(correlation_audit, file.path(table_dir, "sample_correlation_audit.csv"), row.names = FALSE)

all_expression <- do.call(cbind, lapply(names(expression_mats), function(tissue) {
  m <- expression_mats[[tissue]][common_isoforms, , drop = FALSE]
  colnames(m) <- paste(tissue, colnames(m), sep = "__")
  m
}))
log_expression <- log2(all_expression + 1)
row_variance <- apply(log_expression, 1, var, na.rm = TRUE)
top_variable <- names(sort(row_variance, decreasing = TRUE))[seq_len(min(1500, length(row_variance)))]
pca <- prcomp(t(log_expression[top_variable, , drop = FALSE]), center = TRUE, scale. = FALSE)
pca_data <- data.frame(
  sample_id = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2]
) |>
  separate(sample_id, into = c("tissue", "sample"), sep = "__", remove = FALSE) |>
  separate(sample, into = c("stage", "replicate"), sep = "_", remove = FALSE) |>
  mutate(stage_label = factor(stage_labels[stage], levels = unname(stage_labels[stage_order])))
write.csv(pca_data, file.path(table_dir, "sample_pca_coordinates.csv"), row.names = FALSE)

variance_explained <- 100 * (pca$sdev^2 / sum(pca$sdev^2))
e155_midbrain_pca <- pca_data |>
  filter(tissue == "Midbrain", stage_label == "E15.5") |>
  mutate(
    direct_label = paste0("Midbrain E15.5 r", replicate),
    label_hjust = 1.05
  )
p_pca <- ggplot(pca_data, aes(PC1, PC2, color = tissue, shape = stage_label)) +
  geom_point(size = 2.5, alpha = 0.9) +
  geom_text(
    data = e155_midbrain_pca,
    aes(PC1, PC2, label = direct_label, hjust = label_hjust),
    inherit.aes = FALSE,
    color = "#4D4D4D",
    size = 2.4,
    vjust = -0.55,
    show.legend = FALSE
  ) +
  scale_color_manual(values = tissue_colors, name = "Region") +
  scale_shape_manual(values = c(16, 17, 15, 18, 8, 3, 7, 4), name = "Stage") +
  labs(
    title = "Archive-level sample-expression PCA",
    subtitle = sprintf(
      "PCA of the 1,500 most variable isoforms among %s common background isoforms",
      format(length(common_isoforms), big.mark = ",")
    ),
    x = sprintf("PC1 (%.1f%%)", variance_explained[1]),
    y = sprintf("PC2 (%.1f%%)", variance_explained[2]),
    shape = "Stage"
  ) +
  theme_paper(9.5) +
  theme(legend.position = "right")
save_figure(p_pca, "figureS1_sample_pca", 7.2, 5.4)

# Machine-readable summary and software record -------------------------------
summary_values <- data.frame(
  metric = c(
    "Total libraries", "Stages per region", "Pairwise contrasts per region",
    "Analyzable gene universe", "DTU gene union", "Shared DTU core", "Midbrain-only DTU genes",
    "Forebrain DTU genes", "Hindbrain DTU genes", "Midbrain DTU genes"
  ),
  value = c(
    sum(dataset_summary$Libraries), 8, choose(8, 2),
    length(analysis_universe), length(gene_union), length(gene_core), length(midbrain_only),
    length(gene_sets$Forebrain), length(gene_sets$Hindbrain), length(gene_sets$Midbrain)
  )
)
write.csv(summary_values, file.path(table_dir, "manuscript_numbers.csv"), row.names = FALSE)

capture.output(sessionInfo(), file = file.path(data_dir, "sessionInfo.txt"))
cat("Analysis complete. Outputs written to", normalizePath(paper_dir), "\n")
