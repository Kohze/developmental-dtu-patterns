#!/usr/bin/env Rscript

# Test whether the E15.5 episode set is dominated by transcript-end events or
# a systematic shift toward shorter annotated transcripts. These patterns
# would be compatible with sample-wide RNA-processing or quantification bias
# and therefore belong in the counter-audit of the biological interpretation.

suppressPackageStartupMessages({
  library(Biostrings)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else "."
paper_dir <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/")
input_dir <- if (length(args)) {
  normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
} else {
  normalizePath(file.path(paper_dir, "..", ".."), winslash = "/", mustWork = TRUE)
}

load(file.path(input_dir, "isoform_final.RData"))

table_dir <- file.path(paper_dir, "tables")
figure_dir <- file.path(paper_dir, "figures")
episode_path <- file.path(paper_dir, "data", "transient_regional_isoform_episodes.csv")
episodes <- read.csv(episode_path, stringsAsFactors = FALSE, check.names = FALSE)

features <- as.data.frame(combinedMitbrain$isoformFeatures)
is_primary_dtu <- !is.na(features$isoform_switch_q_value) &
  features$isoform_switch_q_value < 0.05 &
  abs(features$dIF) >= 0.10
midbrain_dtu_ids <- unique(features$isoform_id[is_primary_dtu])

event_names <- c("ES", "MEE", "MES", "IR", "A5", "A3", "ATSS", "ATTS")
event_labels <- c(
  ES = "Exon skipping",
  MEE = "Mutually exclusive exon",
  MES = "Multiple exon skipping",
  IR = "Intron retention",
  A5 = "Alternative 5' splice site",
  A3 = "Alternative 3' splice site",
  ATSS = "Alternative transcript start",
  ATTS = "Alternative transcript end"
)
event_data <- as.data.frame(combinedMitbrain$AlternativeSplicingAnalysis)[
  , c("isoform_id", event_names)
]
for (event in event_names) {
  event_data[[event]] <- !is.na(event_data[[event]]) & event_data[[event]] != 0
}

sequence_lengths <- data.frame(
  isoform_id = names(combinedMitbrain$ntSequence),
  transcript_length_nt = as.integer(width(combinedMitbrain$ntSequence)),
  stringsAsFactors = FALSE
)
isoform_map <- features |>
  dplyr::distinct(isoform_id, gene_id)
event_data <- event_data |>
  left_join(isoform_map, by = "isoform_id")

e155 <- episodes |>
  filter(start_stage == 15.5, end_stage == 15.5, replicate_consistent) |>
  distinct(isoform_id, gene_id, gene_name, direction, .keep_all = TRUE)
other_episodes <- episodes |>
  filter(!(start_stage == 15.5 & end_stage == 15.5), replicate_consistent) |>
  distinct(isoform_id)

cohort_ids <- list(
  "All midbrain DTU isoforms" = midbrain_dtu_ids,
  "E15.5 episode isoforms" = e155$isoform_id,
  "Other episode isoforms" = other_episodes$isoform_id
)

event_summary <- bind_rows(lapply(names(cohort_ids), function(cohort) {
  ids <- cohort_ids[[cohort]]
  value <- event_data[event_data$isoform_id %in% ids, , drop = FALSE]
  bind_rows(lapply(event_names, function(event) {
    data.frame(
      cohort = cohort,
      event = event,
      event_label = unname(event_labels[[event]]),
      positive = sum(value[[event]], na.rm = TRUE),
      total = nrow(value),
      prevalence = mean(value[[event]], na.rm = TRUE)
    )
  }))
}))

all_dtu_events <- event_data[event_data$isoform_id %in% midbrain_dtu_ids, ]
e155_events <- event_data[event_data$isoform_id %in% e155$isoform_id, ]
comparison_events <- all_dtu_events[
  !all_dtu_events$isoform_id %in% e155_events$isoform_id,
]
event_tests <- bind_rows(lapply(event_names, function(event) {
  e_pos <- sum(e155_events[[event]])
  b_pos <- sum(comparison_events[[event]])
  cluster_input <- bind_rows(
    e155_events |>
      transmute(gene_id, cohort = "episode", positive = as.integer(.data[[event]])),
    comparison_events |>
      transmute(gene_id, cohort = "comparison", positive = as.integer(.data[[event]]))
  ) |>
    group_by(gene_id, cohort) |>
    summarise(positive = sum(positive), total = n(), .groups = "drop") |>
    pivot_wider(
      names_from = cohort,
      values_from = c(positive, total),
      values_fill = 0
    )
  set.seed(20260730 + match(event, event_names))
  bootstrap_difference <- replicate(5000, {
    sampled <- cluster_input[
      sample.int(nrow(cluster_input), replace = TRUE),
      ,
      drop = FALSE
    ]
    sum(sampled$positive_episode) / sum(sampled$total_episode) -
      sum(sampled$positive_comparison) / sum(sampled$total_comparison)
  })
  data.frame(
    event = event,
    event_label = unname(event_labels[[event]]),
    e155_positive = e_pos,
    e155_total = nrow(e155_events),
    e155_prevalence = e_pos / nrow(e155_events),
    comparison_positive = b_pos,
    comparison_total = nrow(comparison_events),
    comparison_prevalence = b_pos / nrow(comparison_events),
    prevalence_difference = e_pos / nrow(e155_events) - b_pos / nrow(comparison_events),
    cluster_bootstrap_ci_low = unname(quantile(bootstrap_difference, 0.025)),
    cluster_bootstrap_ci_high = unname(quantile(bootstrap_difference, 0.975))
  )
}))

length_data <- sequence_lengths |>
  mutate(
    cohort = case_when(
      isoform_id %in% e155$isoform_id ~ "E15.5 episode isoforms",
      isoform_id %in% midbrain_dtu_ids ~ "Other midbrain DTU isoforms",
      TRUE ~ NA_character_
    )
  ) |>
  filter(!is.na(cohort))

same_gene_non_episode <- isoform_map |>
  filter(
    gene_id %in% e155$gene_id,
    isoform_id %in% midbrain_dtu_ids,
    !isoform_id %in% e155$isoform_id
  ) |>
  left_join(sequence_lengths, by = "isoform_id") |>
  filter(!is.na(transcript_length_nt))

direction_lengths <- e155 |>
  dplyr::select(isoform_id, gene_id, gene_name, direction) |>
  left_join(sequence_lengths, by = "isoform_id") |>
  filter(!is.na(transcript_length_nt))

paired_length <- direction_lengths |>
  group_by(gene_id, gene_name) |>
  summarise(
    n_higher = sum(direction == "higher"),
    n_lower = sum(direction == "lower"),
    mean_higher_length_nt = mean(transcript_length_nt[direction == "higher"]),
    mean_lower_length_nt = mean(transcript_length_nt[direction == "lower"]),
    higher_minus_lower_nt = mean_higher_length_nt - mean_lower_length_nt,
    .groups = "drop"
  ) |>
  filter(n_higher > 0, n_lower > 0, is.finite(higher_minus_lower_nt))

length_summary <- bind_rows(
  length_data |>
    group_by(cohort) |>
    summarise(
      metric = "cohort transcript length",
      n = n(),
      median = median(transcript_length_nt),
      q25 = quantile(transcript_length_nt, 0.25),
      q75 = quantile(transcript_length_nt, 0.75),
      .groups = "drop"
    ),
  direction_lengths |>
    group_by(direction) |>
    summarise(
      metric = "direction-specific transcript length",
      n = n(),
      median = stats::median(transcript_length_nt),
      q25 = stats::quantile(transcript_length_nt, 0.25),
      q75 = stats::quantile(transcript_length_nt, 0.75),
      .groups = "drop"
    ) |>
    mutate(cohort = paste("E15.5", direction)) |>
    dplyr::select(cohort, metric, n, median, q25, q75),
  same_gene_non_episode |>
    summarise(
      cohort = "Non-episode DTU isoforms in episode genes",
      metric = "same-gene comparison transcript length",
      n = n(),
      median = stats::median(transcript_length_nt),
      q25 = stats::quantile(transcript_length_nt, 0.25),
      q75 = stats::quantile(transcript_length_nt, 0.75)
    ),
  data.frame(
    cohort = "Reciprocal genes",
    metric = "higher-minus-lower transcript length",
    n = nrow(paired_length),
    median = median(paired_length$higher_minus_lower_nt),
    q25 = quantile(paired_length$higher_minus_lower_nt, 0.25),
    q75 = quantile(paired_length$higher_minus_lower_nt, 0.75)
  )
)

direction_test <- wilcox.test(
  transcript_length_nt ~ direction,
  data = direction_lengths,
  exact = FALSE
)
paired_sign_test <- binom.test(
  sum(paired_length$higher_minus_lower_nt < 0),
  sum(paired_length$higher_minus_lower_nt != 0),
  p = 0.5
)
same_gene_length_test <- wilcox.test(
  direction_lengths$transcript_length_nt,
  same_gene_non_episode$transcript_length_nt,
  exact = FALSE
)
test_summary <- data.frame(
  test = c(
    "E15.5 higher-versus-lower isoform length (Wilcoxon rank-sum)",
    "Reciprocal genes with higher isoforms shorter than lower isoforms (sign test)",
    "E15.5 episode versus non-episode DTU isoform length within episode genes (Wilcoxon rank-sum)"
  ),
  estimate = c(
    median(direction_lengths$transcript_length_nt[direction_lengths$direction == "higher"]) -
      median(direction_lengths$transcript_length_nt[direction_lengths$direction == "lower"]),
    unname(paired_sign_test$estimate),
    median(direction_lengths$transcript_length_nt) -
      median(same_gene_non_episode$transcript_length_nt)
  ),
  p_value = c(
    direction_test$p.value,
    paired_sign_test$p.value,
    same_gene_length_test$p.value
  ),
  n = c(
    nrow(direction_lengths),
    nrow(paired_length),
    nrow(direction_lengths) + nrow(same_gene_non_episode)
  )
)

write.csv(
  event_tests,
  file.path(table_dir, "e155_episode_event_bias.csv"),
  row.names = FALSE
)
write.csv(
  length_summary,
  file.path(table_dir, "e155_episode_length_summary.csv"),
  row.names = FALSE
)
write.csv(
  paired_length,
  file.path(table_dir, "e155_episode_length_by_gene.csv"),
  row.names = FALSE
)
write.csv(
  test_summary,
  file.path(table_dir, "e155_episode_artifact_tests.csv"),
  row.names = FALSE
)

plot_events <- event_summary |>
  filter(cohort != "Other episode isoforms") |>
  mutate(
    event_label = factor(event_label, levels = rev(unname(event_labels))),
    cohort = factor(
      cohort,
      levels = c("All midbrain DTU isoforms", "E15.5 episode isoforms")
    )
  )

p_event <- ggplot(plot_events, aes(prevalence, event_label, colour = cohort, shape = cohort)) +
  geom_point(size = 2.8, position = position_dodge(width = 0.35)) +
  scale_colour_manual(values = c("#7A8793", "#D55E00")) +
  scale_shape_manual(values = c(16, 17)) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 0.55)) +
  labs(
    title = "a   Exon skipping is enriched in\nE15.5 episodes",
    x = "Fraction carrying annotation",
    y = NULL,
    colour = NULL,
    shape = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

p_length <- ggplot(
  direction_lengths,
  aes(direction, transcript_length_nt, fill = direction)
) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.75, colour = NA) +
  geom_boxplot(width = 0.16, outlier.shape = NA, fill = "white", linewidth = 0.35) +
  scale_fill_manual(values = c(higher = "#0072B2", lower = "#D55E00")) +
  scale_y_log10(labels = scales::label_number(big.mark = ",")) +
  labs(
    title = "b   Similar transcript lengths\nby episode direction",
    x = "Midbrain episode direction",
    y = "Annotated transcript length (nt)",
    fill = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

figure <- p_event + p_length + plot_layout(widths = c(1.15, 0.85))
ggsave(
  file.path(figure_dir, "figureS4_e155_artifact_audit.pdf"),
  figure,
  width = 9.2,
  height = 4.8,
  device = cairo_pdf
)
ggsave(
  file.path(figure_dir, "figureS4_e155_artifact_audit.png"),
  figure,
  width = 9.2,
  height = 4.8,
  dpi = 300,
  bg = "white"
)

print(event_tests)
print(length_summary)
print(test_summary)
