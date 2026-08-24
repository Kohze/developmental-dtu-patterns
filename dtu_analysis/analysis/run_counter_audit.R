#!/usr/bin/env Rscript

# Adversarial robustness audit for the DTU manuscript.
# Usage:
#   Rscript analysis/run_counter_audit.R [directory containing the two .RData inputs]

options(stringsAsFactors = FALSE, width = 140)
set.seed(20260729)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(stringr)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
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

table_dir <- file.path(paper_dir, "tables")
data_dir <- file.path(paper_dir, "data")
figure_dir <- file.path(paper_dir, "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

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

stage_order <- c("10.5", "11.5", "12.5", "13.5", "14.5", "15.5", "16.5", "0")
stage_labels <- c(
  "10.5" = "E10.5", "11.5" = "E11.5", "12.5" = "E12.5",
  "13.5" = "E13.5", "14.5" = "E14.5", "15.5" = "E15.5",
  "16.5" = "E16.5", "0" = "P0"
)
normalise_condition <- function(x) sub("^X", "", as.character(x))

object_integrity <- bind_rows(lapply(names(objects), function(region) {
  analysed <- objects[[region]]
  background <- background_objects[[region]]
  data.frame(
    region = region,
    archived_package_version = background$runInfo$IsoformSwitchAnalyzeR$version,
    libraries = nrow(background$designMatrix),
    stages = length(unique(background$designMatrix$condition)),
    background_isoforms = nrow(background$isoformRepExpression),
    analysed_isoforms = nrow(analysed$isoformRepExpression),
    analysed_feature_rows = nrow(analysed$isoformFeatures),
    distinct_contrasts = nrow(unique(data.frame(
      condition_1 = normalise_condition(analysed$isoformFeatures$condition_1),
      condition_2 = normalise_condition(analysed$isoformFeatures$condition_2)
    )))
  )
}))
write.csv(object_integrity, file.path(table_dir, "audit_object_integrity.csv"), row.names = FALSE)

qvalue_family <- bind_rows(lapply(names(objects), function(region) {
  analysed <- as.data.frame(objects[[region]]$isoformFeatures)
  background <- as.data.frame(background_objects[[region]]$isoformFeatures)
  input_counts <- background |>
    dplyr::count(condition_1, condition_2, name = "input_isoforms")
  retained_counts <- analysed |>
    dplyr::count(condition_1, condition_2, name = "retained_rows")
  contrast_counts <- left_join(
    input_counts,
    retained_counts,
    by = c("condition_1", "condition_2")
  )
  data.frame(
    region = region,
    contrasts = nrow(contrast_counts),
    input_isoforms_per_contrast_min = min(contrast_counts$input_isoforms),
    input_isoforms_per_contrast_max = max(contrast_counts$input_isoforms),
    retained_rows_per_contrast_min = min(contrast_counts$retained_rows),
    retained_rows_per_contrast_max = max(contrast_counts$retained_rows),
    adjustment_scope = paste(
      "One DEXSeq model and Benjamini-Hochberg-adjusted isoform p-value",
      "family per region-specific stage-pair contrast"
    ),
    across_contrast_control = FALSE,
    across_region_control = FALSE,
    exact_recomputation_from_reduced_object = FALSE
  )
}))
write.csv(
  qvalue_family,
  file.path(table_dir, "audit_qvalue_family.csv"),
  row.names = FALSE
)

background_gene_sets <- lapply(background_objects, function(x) unique(x$isoformFeatures$gene_id))
analysis_universe <- Reduce(union, background_gene_sets)
universe_intersection <- Reduce(intersect, background_gene_sets)

get_gene_set <- function(x, q_cutoff = 0.05, dif_cutoff = 0.10) {
  d <- as.data.frame(x$isoformFeatures)
  keep <- !is.na(d$isoform_switch_q_value) &
    d$isoform_switch_q_value < q_cutoff &
    abs(d$dIF) >= dif_cutoff
  unique(d$gene_id[keep])
}

primary_sets <- lapply(objects, get_gene_set)
primary_union <- Reduce(union, primary_sets)
primary_core <- Reduce(intersect, primary_sets)
primary_midbrain_only <- setdiff(
  primary_sets$Midbrain,
  union(primary_sets$Forebrain, primary_sets$Hindbrain)
)

pairwise_overlap <- bind_rows(lapply(combn(names(primary_sets), 2, simplify = FALSE), function(pair) {
  a <- primary_sets[[pair[1]]]
  b <- primary_sets[[pair[2]]]
  observed <- length(intersect(a, b))
  expected <- length(a) * length(b) / length(analysis_universe)
  p_upper <- phyper(
    observed - 1,
    length(a),
    length(analysis_universe) - length(a),
    length(b),
    lower.tail = FALSE
  )
  data.frame(
    comparison = paste(pair, collapse = " vs "),
    set_1 = length(a),
    set_2 = length(b),
    observed_overlap = observed,
    expected_overlap_independence = expected,
    fold_over_independence = observed / expected,
    jaccard = observed / length(union(a, b)),
    overlap_coefficient = observed / min(length(a), length(b)),
    hypergeometric_p = p_upper
  )
}))
triple_expected <- prod(vapply(primary_sets, length, integer(1))) / length(analysis_universe)^2
overlap_summary <- bind_rows(
  pairwise_overlap,
  data.frame(
    comparison = "All three regions",
    set_1 = NA_integer_,
    set_2 = NA_integer_,
    observed_overlap = length(primary_core),
    expected_overlap_independence = triple_expected,
    fold_over_independence = length(primary_core) / triple_expected,
    jaccard = length(primary_core) / length(primary_union),
    overlap_coefficient = length(primary_core) / min(vapply(primary_sets, length, integer(1))),
    hypergeometric_p = NA_real_
  )
)
write.csv(overlap_summary, file.path(table_dir, "audit_overlap_statistics.csv"), row.names = FALSE)

threshold_grid <- expand.grid(
  q_cutoff = c(0.05, 0.01, 0.05 / 28),
  dIF_cutoff = c(0.10, 0.20),
  KEEP.OUT.ATTRS = FALSE
) |>
  arrange(desc(q_cutoff), dIF_cutoff)

threshold_sensitivity <- bind_rows(lapply(seq_len(nrow(threshold_grid)), function(i) {
  q_cutoff <- threshold_grid$q_cutoff[i]
  dif_cutoff <- threshold_grid$dIF_cutoff[i]
  sets <- lapply(objects, get_gene_set, q_cutoff = q_cutoff, dif_cutoff = dif_cutoff)
  union_set <- Reduce(union, sets)
  core_set <- Reduce(intersect, sets)
  midbrain_only <- setdiff(sets$Midbrain, union(sets$Forebrain, sets$Hindbrain))
  data.frame(
    q_cutoff = q_cutoff,
    dIF_cutoff = dif_cutoff,
    forebrain_genes = length(sets$Forebrain),
    hindbrain_genes = length(sets$Hindbrain),
    midbrain_genes = length(sets$Midbrain),
    union_genes = length(union_set),
    shared_core = length(core_set),
    shared_fraction_union = length(core_set) / length(union_set),
    midbrain_only = length(midbrain_only)
  )
}))
write.csv(
  threshold_sensitivity,
  file.path(table_dir, "audit_threshold_sensitivity.csv"),
  row.names = FALSE
)

count_contrast_genes <- function(x, q_cutoff, dif_cutoff) {
  d <- as.data.frame(x$isoformFeatures)
  d$c1 <- normalise_condition(d$condition_1)
  d$c2 <- normalise_condition(d$condition_2)
  d <- d[
    !is.na(d$isoform_switch_q_value) &
      d$isoform_switch_q_value < q_cutoff &
      abs(d$dIF) >= dif_cutoff,
    ,
    drop = FALSE
  ]
  d$i1 <- match(d$c1, stage_order)
  d$i2 <- match(d$c2, stage_order)
  d$early_i <- pmin(d$i1, d$i2)
  d$late_i <- pmax(d$i1, d$i2)
  d |>
    distinct(gene_id, early_i, late_i) |>
    dplyr::count(early_i, late_i, name = "genes")
}

e155_sensitivity <- bind_rows(lapply(seq_len(nrow(threshold_grid)), function(i) {
  q_cutoff <- threshold_grid$q_cutoff[i]
  dif_cutoff <- threshold_grid$dIF_cutoff[i]
  all_adjacent <- bind_rows(lapply(names(objects), function(region) {
    count_contrast_genes(objects[[region]], q_cutoff, dif_cutoff) |>
      filter(late_i - early_i == 1) |>
      mutate(region = region)
  }))
  midbrain <- all_adjacent |> filter(region == "Midbrain")
  non_flanking <- midbrain |> filter(!(early_i %in% c(5, 6)))
  data.frame(
    q_cutoff = q_cutoff,
    dIF_cutoff = dif_cutoff,
    midbrain_E14_5_E15_5 = midbrain$genes[midbrain$early_i == 5],
    midbrain_E15_5_E16_5 = midbrain$genes[midbrain$early_i == 6],
    median_other_midbrain_adjacent = median(non_flanking$genes),
    E14_5_E15_5_fold_vs_other_median =
      midbrain$genes[midbrain$early_i == 5] / median(non_flanking$genes),
    E15_5_E16_5_fold_vs_other_median =
      midbrain$genes[midbrain$early_i == 6] / median(non_flanking$genes)
  )
}))
write.csv(e155_sensitivity, file.path(table_dir, "audit_e155_sensitivity.csv"), row.names = FALSE)

criterion_labels <- with(
  threshold_sensitivity,
  sprintf(
    "q < %s\n|dIF| >= %.2f",
    ifelse(
      q_cutoff < 0.002,
      "0.05/28",
      ifelse(q_cutoff < 0.02, "0.01", "0.05")
    ),
    dIF_cutoff
  )
)
threshold_plot_data <- threshold_sensitivity |>
  mutate(criterion = factor(criterion_labels, levels = criterion_labels)) |>
  dplyr::select(criterion, shared_core, midbrain_only) |>
  pivot_longer(
    cols = c(shared_core, midbrain_only),
    names_to = "set",
    values_to = "genes"
  ) |>
  mutate(set = recode(
    set,
    shared_core = "Three-region core",
    midbrain_only = "Midbrain-only detection"
  ))

e155_plot_data <- e155_sensitivity |>
  mutate(criterion = factor(criterion_labels, levels = criterion_labels)) |>
  dplyr::select(criterion, midbrain_E14_5_E15_5, midbrain_E15_5_E16_5) |>
  pivot_longer(
    cols = c(midbrain_E14_5_E15_5, midbrain_E15_5_E16_5),
    names_to = "contrast",
    values_to = "genes"
  ) |>
  mutate(contrast = recode(
    contrast,
    midbrain_E14_5_E15_5 = "E14.5-E15.5",
    midbrain_E15_5_E16_5 = "E15.5-E16.5"
  ))

audit_theme <- theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.position = "top",
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1)
  )
robustness_a <- ggplot(
  threshold_plot_data,
  aes(criterion, genes, color = set, shape = set, group = set)
) +
  geom_point(size = 2.6, position = position_dodge(width = 0.28)) +
  geom_text(
    aes(label = genes),
    size = 3.0,
    vjust = -0.8,
    position = position_dodge(width = 0.28),
    show.legend = FALSE
  ) +
  scale_color_manual(values = c(
    "Three-region core" = "#59636F",
    "Midbrain-only detection" = "#D55E00"
  )) +
  scale_shape_manual(values = c(
    "Three-region core" = 16,
    "Midbrain-only detection" = 17
  )) +
  labs(
    title = "A  Regional structure",
    x = "Detection criterion",
    y = "Genes"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.14))) +
  audit_theme
robustness_b <- ggplot(
  e155_plot_data,
  aes(criterion, genes, color = contrast, shape = contrast, group = contrast)
) +
  geom_point(size = 2.6, position = position_dodge(width = 0.50)) +
  geom_text(
    aes(label = genes),
    size = 3.0,
    vjust = -0.8,
    position = position_dodge(width = 0.50),
    show.legend = FALSE
  ) +
  scale_color_manual(values = c(
    "E14.5-E15.5" = "#0072B2",
    "E15.5-E16.5" = "#D55E00"
  )) +
  scale_shape_manual(values = c(
    "E14.5-E15.5" = 16,
    "E15.5-E16.5" = 17
  )) +
  labs(
    title = "B  Midbrain contrasts flanking E15.5",
    x = "Detection criterion",
    y = "DTU genes"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.14))) +
  audit_theme
robustness_figure <- robustness_a + robustness_b +
  plot_annotation(
    title = "Central DTU patterns persist under stricter thresholds",
    subtitle = "The 0.05/28 criterion is a sensitivity threshold, not joint FDR control"
  )
ggsave(
  file.path(figure_dir, "figureS2_threshold_robustness.pdf"),
  robustness_figure,
  width = 8.3,
  height = 5.2,
  units = "in",
  device = cairo_pdf
)
ggsave(
  file.path(figure_dir, "figureS2_threshold_robustness.png"),
  robustness_figure,
  width = 8.3,
  height = 5.2,
  units = "in",
  dpi = 300,
  bg = "white"
)

expression_matrix <- function(x, field = "isoformRepExpression") {
  d <- as.data.frame(x[[field]])
  rownames(d) <- d$isoform_id
  as.matrix(d[, setdiff(names(d), "isoform_id"), drop = FALSE])
}
background_expression <- lapply(background_objects, expression_matrix)
common_background_isoforms <- Reduce(intersect, lapply(background_expression, rownames))

within_region_correlations <- bind_rows(lapply(names(background_expression), function(region) {
  m <- background_expression[[region]][common_background_isoforms, , drop = FALSE]
  bind_rows(lapply(stage_order, function(stage) {
    data.frame(
      comparison = "Within-region biological replicates",
      region_pair = region,
      stage = unname(stage_labels[stage]),
      replicate = "1 versus 2",
      pearson_log1p = cor(
        log1p(m[, paste0(stage, "_1")]),
        log1p(m[, paste0(stage, "_2")]),
        use = "pairwise.complete.obs"
      )
    )
  }))
}))
same_stage_correlations <- bind_rows(lapply(
  list(c("Forebrain", "Hindbrain"), c("Forebrain", "Midbrain"), c("Hindbrain", "Midbrain")),
  function(pair) {
    a <- background_expression[[pair[1]]][common_background_isoforms, , drop = FALSE]
    b <- background_expression[[pair[2]]][common_background_isoforms, , drop = FALSE]
    bind_rows(lapply(stage_order, function(stage) {
      bind_rows(lapply(1:2, function(rep_a) {
        bind_rows(lapply(1:2, function(rep_b) {
          sample_a <- paste0(stage, "_", rep_a)
          sample_b <- paste0(stage, "_", rep_b)
          data.frame(
            comparison = "Same-stage cross-region replicate combinations",
            region_pair = paste(pair, collapse = " vs "),
            stage = unname(stage_labels[stage]),
            replicate = paste(rep_a, "versus", rep_b),
            pearson_log1p = cor(
              log1p(a[, sample_a]),
              log1p(b[, sample_b]),
              use = "pairwise.complete.obs"
            )
          )
        }))
      }))
    }))
  }
))
background_correlation_audit <- bind_rows(
  within_region_correlations,
  same_stage_correlations
)
write.csv(
  background_correlation_audit,
  file.path(table_dir, "sample_correlation_audit.csv"),
  row.names = FALSE
)

all_background_expression <- do.call(cbind, lapply(names(background_expression), function(region) {
  m <- background_expression[[region]][common_background_isoforms, , drop = FALSE]
  colnames(m) <- paste(region, colnames(m), sep = "__")
  m
}))
log_background_expression <- log2(all_background_expression + 1)
row_variance <- apply(log_background_expression, 1, var, na.rm = TRUE)
top_variable <- names(sort(row_variance, decreasing = TRUE))[
  seq_len(min(1500, length(row_variance)))
]
pca <- prcomp(
  t(log_background_expression[top_variable, , drop = FALSE]),
  center = TRUE,
  scale. = FALSE
)
pca_data <- data.frame(
  sample_id = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2]
) |>
  separate(sample_id, into = c("region", "sample"), sep = "__", remove = FALSE) |>
  separate(sample, into = c("stage", "replicate"), sep = "_", remove = FALSE) |>
  mutate(stage_label = unname(stage_labels[stage]))
variance_explained <- 100 * pca$sdev^2 / sum(pca$sdev^2)
pca_data$PC1_variance_percent <- variance_explained[1]
pca_data$PC2_variance_percent <- variance_explained[2]
write.csv(pca_data, file.path(table_dir, "sample_pca_coordinates.csv"), row.names = FALSE)

midbrain_features <- as.data.frame(objects$Midbrain$isoformFeatures)
midbrain_features$c1 <- normalise_condition(midbrain_features$condition_1)
midbrain_features$c2 <- normalise_condition(midbrain_features$condition_2)
midbrain_features$significant <- !is.na(midbrain_features$isoform_switch_q_value) &
  midbrain_features$isoform_switch_q_value < 0.05 &
  abs(midbrain_features$dIF) >= 0.10

has_significant_contrast <- function(gene, a, b) {
  d <- midbrain_features[
    midbrain_features$gene_id == gene &
      midbrain_features$significant &
      ((midbrain_features$c1 == a & midbrain_features$c2 == b) |
         (midbrain_features$c1 == b & midbrain_features$c2 == a)),
    ,
    drop = FALSE
  ]
  nrow(d) > 0
}

midbrain_if <- expression_matrix(objects$Midbrain, "isoformRepIF")
midbrain_expression <- expression_matrix(objects$Midbrain, "isoformRepExpression")
candidate_rows <- list()
for (gene in primary_midbrain_only) {
  if (!has_significant_contrast(gene, "14.5", "15.5") ||
      !has_significant_contrast(gene, "15.5", "16.5")) {
    next
  }
  ids <- unique(midbrain_features$isoform_id[midbrain_features$gene_id == gene])
  ids <- intersect(ids, intersect(rownames(midbrain_if), rownames(midbrain_expression)))
  if (!length(ids)) next
  gene_tpm_e155 <- sum(midbrain_expression[ids, c("15.5_1", "15.5_2"), drop = FALSE]) / 2
  for (id in ids) {
    values <- as.numeric(midbrain_if[
      id,
      c("14.5_1", "14.5_2", "15.5_1", "15.5_2", "16.5_1", "16.5_2")
    ])
    if (anyNA(values)) next
    e155 <- mean(values[3:4])
    flank_mean <- mean(values[c(1, 2, 5, 6)])
    candidate_rows[[length(candidate_rows) + 1]] <- data.frame(
      gene = gene,
      isoform_id = id,
      transient_IF_score = abs(e155 - flank_mean),
      E15_5_mean_IF = e155,
      flanking_mean_IF = flank_mean,
      E15_5_replicate_gap = abs(values[3] - values[4]),
      gene_TPM_E15_5 = gene_tpm_e155
    )
  }
}
candidate_isoforms <- bind_rows(candidate_rows)
candidate_ranking <- candidate_isoforms |>
  filter(gene_TPM_E15_5 >= 5, E15_5_replicate_gap <= 0.15) |>
  group_by(gene) |>
  summarise(
    max_transient_IF_score = max(transient_IF_score),
    gene_TPM_E15_5 = dplyr::first(gene_TPM_E15_5),
    best_isoform = isoform_id[which.max(transient_IF_score)],
    best_isoform_replicate_gap = E15_5_replicate_gap[which.max(transient_IF_score)],
    .groups = "drop"
  ) |>
  arrange(desc(max_transient_IF_score), best_isoform_replicate_gap)
write.csv(candidate_ranking, file.path(table_dir, "audit_candidate_ranking.csv"), row.names = FALSE)

# Candidate accession and archived-annotation audit. Current RefSeq versions
# are cross-references only: the archived objects retain unversioned accessions,
# and the original annotation release is not available.
candidate_current_refseq <- data.frame(
  gene = c("Ppp2r3a", "Ppp2r3a", "Rtn2", "Rtn2"),
  archived_isoform_id = c("NM_001161362", "NM_172144", "NM_001025364", "NM_013648"),
  ncbi_gene_id = c(235542L, 235542L, 20167L, 20167L),
  current_refseq_version = c("NM_001161362.3", "NM_172144.3", "NM_001025364.3", "NM_013648.6"),
  current_protein_accession = c("NP_001154834.1", "NP_742156.2", "NP_001020535.1", "NP_038676.1"),
  current_transcript_label = c("variant 1 / isoform 1", "variant 2 / isoform 2",
                               "variant C / isoform C", "variant B / isoform B"),
  current_refseq_status = "VALIDATED",
  current_check_date = "2026-07-30",
  current_source_url = c(
    "https://www.ncbi.nlm.nih.gov/gene/235542",
    "https://www.ncbi.nlm.nih.gov/gene/235542",
    "https://www.ncbi.nlm.nih.gov/gene/20167",
    "https://www.ncbi.nlm.nih.gov/gene/20167"
  )
)

candidate_exons <- as.data.frame(objects$Midbrain$exons)
candidate_orfs <- as.data.frame(objects$Midbrain$orfAnalysis)
candidate_splicing <- as.data.frame(objects$Midbrain$AlternativeSplicingAnalysis)
candidate_feature_map <- midbrain_features |>
  distinct(gene_id, isoform_id)

candidate_identifier_rows <- lapply(seq_len(nrow(candidate_current_refseq)), function(i) {
  row <- candidate_current_refseq[i, , drop = FALSE]
  id <- row$archived_isoform_id
  gene <- row$gene
  exon_rows <- candidate_exons[candidate_exons$isoform_id == id, , drop = FALSE]
  orf_row <- candidate_orfs[candidate_orfs$isoform_id == id, , drop = FALSE]
  splicing_row <- candidate_splicing[candidate_splicing$isoform_id == id, , drop = FALSE]
  nt_length <- nchar(as.character(objects$Midbrain$ntSequence[id]))
  aa_length <- nchar(as.character(objects$Midbrain$aaSequence[id]))
  data.frame(
    row,
    archived_gene_mapping = unique(exon_rows$gene_id),
    archived_isoforms_for_gene = n_distinct(
      candidate_feature_map$isoform_id[candidate_feature_map$gene_id == gene]
    ),
    archived_exon_count = nrow(exon_rows),
    archived_spliced_length_nt = sum(exon_rows$width),
    archived_nt_sequence_length = nt_length,
    archived_orf_length_nt = orf_row$orfTransciptLength,
    archived_protein_length_aa = aa_length,
    archived_PTC = orf_row$PTC,
    archived_ES_flag = splicing_row$ES,
    archived_ATSS_flag = splicing_row$ATSS,
    base_accession_matches_current = identical(
      id, sub("\\.[0-9]+$", "", row$current_refseq_version)
    ),
    archived_version_status = paste(
      "Unresolved: archived identifier omits the version suffix and",
      "the original annotation release is unavailable"
    )
  )
})
candidate_identifier_audit <- bind_rows(candidate_identifier_rows)
write.csv(
  candidate_identifier_audit,
  file.path(table_dir, "candidate_identifier_audit.csv"),
  row.names = FALSE
)

# Functional-enrichment counter-audits --------------------------------------
# Genes with more detectable isoforms and higher expression have more
# opportunities to enter a DTU set. Test whether the leading shared-core GO
# terms retain an association after adjustment for those two archive-derived
# detection-opportunity covariates.
gene_detection_features <- bind_rows(lapply(names(background_objects), function(region) {
  d <- as.data.frame(background_objects[[region]]$isoformFeatures) |>
    dplyr::select(gene_id, isoform_id, gene_overall_mean) |>
    distinct()
  d |>
    group_by(gene_id) |>
    summarise(
      isoform_count = n_distinct(isoform_id),
      mean_expression = median(gene_overall_mean, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(region = region)
})) |>
  group_by(gene_id) |>
  summarise(
    isoform_count = max(isoform_count, na.rm = TRUE),
    mean_expression = mean(mean_expression, na.rm = TRUE),
    .groups = "drop"
  ) |>
  filter(gene_id %in% analysis_universe) |>
  mutate(
    in_shared_core = as.integer(gene_id %in% primary_core),
    log_isoform_count = log1p(isoform_count),
    log_mean_expression = log1p(pmax(mean_expression, 0))
  )

go_top_path <- file.path(table_dir, "go_bp_enrichment_top.csv")
if (!file.exists(go_top_path)) {
  stop("Run run_analysis.R before the counter-audit; missing ", go_top_path)
}
shared_top_go <- read.csv(go_top_path, check.names = FALSE) |>
  filter(group == "Shared core") |>
  arrange(p.adjust) |>
  slice_head(n = 8)

go_membership <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = analysis_universe,
  keytype = "SYMBOL",
  columns = c("GOALL", "ONTOLOGYALL")
) |>
  filter(ONTOLOGYALL == "BP", !is.na(GOALL)) |>
  distinct(SYMBOL, GOALL)

adjusted_go <- bind_rows(lapply(seq_len(nrow(shared_top_go)), function(i) {
  go_id <- shared_top_go$ID[i]
  model_data <- gene_detection_features |>
    mutate(in_term = as.integer(gene_id %in% go_membership$SYMBOL[go_membership$GOALL == go_id]))
  fit <- glm(
    in_shared_core ~ in_term + log_isoform_count + log_mean_expression,
    data = model_data,
    family = binomial()
  )
  estimate <- summary(fit)$coefficients["in_term", "Estimate"]
  standard_error <- summary(fit)$coefficients["in_term", "Std. Error"]
  p_value <- summary(fit)$coefficients["in_term", "Pr(>|z|)"]
  data.frame(
    ID = go_id,
    Description = shared_top_go$Description[i],
    unadjusted_fold_enrichment = shared_top_go$enrichment[i],
    unadjusted_FDR = shared_top_go$p.adjust[i],
    adjusted_odds_ratio = exp(estimate),
    adjusted_CI_low = exp(estimate - 1.96 * standard_error),
    adjusted_CI_high = exp(estimate + 1.96 * standard_error),
    adjusted_p_value = p_value,
    term_genes_in_universe = sum(model_data$in_term),
    shared_core_genes_in_term = sum(model_data$in_term & model_data$in_shared_core)
  )
}))
adjusted_go$adjusted_FDR <- p.adjust(adjusted_go$adjusted_p_value, method = "BH")
write.csv(
  adjusted_go,
  file.path(table_dir, "audit_go_covariate_adjusted.csv"),
  row.names = FALSE
)

# "Midbrain only" is a threshold-defined detection category. Compare it with
# other DTU genes, rather than only with the full expressed-gene background.
comparative_midbrain_go <- bind_rows(
  lapply(
    list(
    `All DTU genes` = primary_union,
    `All midbrain DTU genes` = primary_sets$Midbrain
    ),
    function(comparison_universe) {
      ego <- enrichGO(
        gene = primary_midbrain_only,
        universe = comparison_universe,
        OrgDb = org.Mm.eg.db,
        keyType = "SYMBOL",
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 1,
        qvalueCutoff = 1,
        minGSSize = 10,
        maxGSSize = 500,
        readable = FALSE
      )
      as.data.frame(ego)
    }
  ),
  .id = "comparison_universe"
)
write.csv(
  comparative_midbrain_go,
  file.path(table_dir, "audit_midbrain_go_comparative.csv"),
  row.names = FALSE
)

# Replace the initial absolute-enrichment panel with a two-part sensitivity
# figure for the shared core. The right panel is a post hoc covariate analysis
# of the eight terms selected from the primary over-representation results.
go_order <- rev(shared_top_go$Description)
plot_go_unadjusted <- shared_top_go |>
  mutate(
    Description = factor(Description, levels = go_order),
    label = str_wrap(as.character(Description), width = 35)
  ) |>
  ggplot(aes(enrichment, Description, size = Count, color = -log10(p.adjust))) +
  geom_point(alpha = 0.9) +
  scale_color_viridis_c(option = "C", end = 0.88, name = expression(-log[10]("FDR"))) +
  scale_size_continuous(range = c(2.5, 7), name = "Genes") +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 35)) +
  labs(
    title = "A  Standard\nover-representation",
    x = "Fold enrichment",
    y = NULL
  ) +
  audit_theme +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

plot_go_adjusted <- adjusted_go |>
  mutate(Description = factor(Description, levels = go_order)) |>
  ggplot(aes(adjusted_odds_ratio, Description, color = adjusted_FDR < 0.05)) +
  geom_vline(xintercept = 1, linetype = 2, color = "#888888") +
  geom_errorbarh(
    aes(xmin = adjusted_CI_low, xmax = adjusted_CI_high),
    height = 0.20,
    linewidth = 0.55
  ) +
  geom_point(size = 2.7) +
  scale_x_log10() +
  scale_color_manual(
    values = c(`TRUE` = "#0072B2", `FALSE` = "#999999"),
    labels = c(`TRUE` = "FDR < 0.05", `FALSE` = "FDR >= 0.05"),
    name = NULL
  ) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 35)) +
  labs(
    title = "B  Post hoc covariate\nsensitivity",
    x = "Adjusted odds ratio (log scale)",
    y = NULL
  ) +
  audit_theme +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "bottom"
  )

figure4_audit <- plot_go_unadjusted + plot_go_adjusted +
  plot_annotation(
    title = "Shared-core GO enrichment and post hoc covariate sensitivity",
    subtitle = "Adjustment covariates: archive-derived isoform count and mean expression"
  )
ggsave(
  file.path(figure_dir, "figure4_go_enrichment.pdf"),
  figure4_audit,
  width = 9.2,
  height = 6.0,
  units = "in",
  device = cairo_pdf
)
ggsave(
  file.path(figure_dir, "figure4_go_enrichment.png"),
  figure4_audit,
  width = 9.2,
  height = 6.0,
  units = "in",
  dpi = 300,
  bg = "white"
)

audit_summary <- data.frame(
  metric = c(
    "Union background genes",
    "Intersection background genes",
    "Common background isoforms used for QC",
    "Primary DTU union",
    "Primary shared core",
    "Primary midbrain-only set",
    "Triple-overlap fold over independence",
    "Smallest comparative midbrain-only GO FDR",
    "Adjusted shared-core GO terms at FDR < 0.05",
    "PCA PC1 variance percent",
    "PCA PC2 variance percent"
  ),
  value = c(
    length(analysis_universe),
    length(universe_intersection),
    length(common_background_isoforms),
    length(primary_union),
    length(primary_core),
    length(primary_midbrain_only),
    length(primary_core) / triple_expected,
    min(comparative_midbrain_go$p.adjust, na.rm = TRUE),
    sum(adjusted_go$adjusted_FDR < 0.05, na.rm = TRUE),
    variance_explained[1],
    variance_explained[2]
  )
)
write.csv(audit_summary, file.path(table_dir, "counter_audit_summary.csv"), row.names = FALSE)

cat("Counter-audit complete.\n")
print(audit_summary, row.names = FALSE)
cat("\nThreshold sensitivity:\n")
print(threshold_sensitivity, row.names = FALSE)
cat("\nE15.5 sensitivity:\n")
print(e155_sensitivity, row.names = FALSE)
cat("\nTop ranked transient candidates:\n")
print(head(candidate_ranking, 10), row.names = FALSE)
