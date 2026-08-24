#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 160)
set.seed(20260730)

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
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
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) {
  stop("Expected at most one argument: the directory containing the two .RData inputs.")
}
input_candidates <- if (length(args)) {
  args[[1]]
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
  Forebrain = switchListForebrain,
  Hindbrain = switchListHindbrain,
  Midbrain = switchListMidbrain
)
regions <- names(objects)
stage_levels <- c("10.5", "11.5", "12.5", "13.5", "14.5", "15.5", "16.5", "0")
stage_codes <- c(
  "10.5" = "E10_5", "11.5" = "E11_5", "12.5" = "E12_5",
  "13.5" = "E13_5", "14.5" = "E14_5", "15.5" = "E15_5",
  "16.5" = "E16_5", "0" = "P0"
)
stage_labels <- c(
  "10.5" = "E10.5", "11.5" = "E11.5", "12.5" = "E12.5",
  "13.5" = "E13.5", "14.5" = "E14.5", "15.5" = "E15.5",
  "16.5" = "E16.5", "0" = "P0"
)

extract_matrix <- function(x, field) {
  value <- x[[field]]
  stopifnot(!anyDuplicated(value$isoform_id))
  answer <- as.matrix(value[, -1, drop = FALSE])
  storage.mode(answer) <- "numeric"
  rownames(answer) <- value$isoform_id
  answer
}

extract_mapping <- function(x) {
  value <- unique(x$isoformFeatures[, c("isoform_id", "gene_id", "gene_name")])
  value[!duplicated(value$isoform_id), , drop = FALSE]
}

count_matrices <- lapply(objects, extract_matrix, field = "isoformCountMatrix")
if_matrices <- lapply(objects, extract_matrix, field = "isoformRepIF")
mapping_tables <- lapply(objects, extract_mapping)
common_isoforms <- Reduce(intersect, lapply(count_matrices, rownames))
mapping <- mapping_tables$Forebrain[
  match(common_isoforms, mapping_tables$Forebrain$isoform_id),
  c("isoform_id", "gene_id", "gene_name")
]
stopifnot(!anyNA(mapping$gene_id))

combine_regions <- function(values) {
  do.call(cbind, lapply(names(values), function(region) {
    answer <- values[[region]][common_isoforms, , drop = FALSE]
    colnames(answer) <- paste(region, colnames(answer), sep = "__")
    answer
  }))
}
counts <- combine_regions(count_matrices)
isoform_fraction <- combine_regions(if_matrices)

parts <- strsplit(colnames(counts), "__", fixed = TRUE)
metadata <- data.frame(
  sample_id = colnames(counts),
  region = vapply(parts, `[[`, character(1), 1),
  stage_replicate = vapply(parts, `[[`, character(1), 2),
  stringsAsFactors = FALSE
)
metadata$stage <- sub("_[12]$", "", metadata$stage_replicate)
metadata$replicate <- sub("^.*_", "", metadata$stage_replicate)
metadata$stage <- factor(metadata$stage, levels = stage_levels)
metadata$region <- factor(metadata$region, levels = regions)
metadata$stage_code <- unname(stage_codes[as.character(metadata$stage)])
group_levels <- unlist(lapply(
  regions,
  function(region) paste(region, unname(stage_codes[stage_levels]), sep = "_")
))
metadata$group <- factor(
  paste(metadata$region, metadata$stage_code, sep = "_"),
  levels = group_levels
)

keep <- filterByExpr(
  counts,
  group = metadata$group,
  min.count = 10,
  min.total.count = 30,
  large.n = 20,
  min.prop = 0.7
)
mapping <- mapping[keep, , drop = FALSE]
multi_isoform <- duplicated(mapping$gene_id) | duplicated(mapping$gene_id, fromLast = TRUE)
mapping <- mapping[multi_isoform, , drop = FALSE]
isoform_fraction <- isoform_fraction[mapping$isoform_id, , drop = FALSE]

design <- model.matrix(~ 0 + group, metadata)
colnames(design) <- levels(metadata$group)
stopifnot(qr(design)$rank == ncol(design))

region_pairs <- combn(regions, 2, simplify = FALSE)
contrast_matrix <- matrix(
  0,
  nrow = ncol(design),
  ncol = length(region_pairs) * length(stage_levels),
  dimnames = list(
    colnames(design),
    character(length(region_pairs) * length(stage_levels))
  )
)
contrast_metadata <- data.frame()
index <- 0L
for (pair in region_pairs) {
  for (stage in stage_levels) {
    index <- index + 1L
    contrast_name <- paste(pair[[1]], "vs", pair[[2]], stage_codes[[stage]], sep = "_")
    colnames(contrast_matrix)[[index]] <- contrast_name
    contrast_matrix[paste(pair[[1]], stage_codes[[stage]], sep = "_"), index] <- 1
    contrast_matrix[paste(pair[[2]], stage_codes[[stage]], sep = "_"), index] <- -1
    contrast_metadata <- rbind(
      contrast_metadata,
      data.frame(
        contrast = contrast_name,
        region_1 = pair[[1]],
        region_2 = pair[[2]],
        stage = stage,
        stringsAsFactors = FALSE
      )
    )
  }
}
stopifnot(all(abs(colSums(contrast_matrix)) < 1e-12))

epsilon <- 0.005
logit_if <- qlogis(pmin(pmax(isoform_fraction, epsilon), 1 - epsilon))
fit <- lmFit(logit_if, design)
fit <- contrasts.fit(fit, contrast_matrix)
fit <- eBayes(fit, trend = TRUE, robust = TRUE)
raw_fit <- contrasts.fit(lmFit(isoform_fraction, design), contrast_matrix)
group_means <- lmFit(isoform_fraction, design)$coefficients

pair_result <- do.call(rbind, lapply(seq_len(ncol(contrast_matrix)), function(i) {
  data.frame(
    isoform_id = rownames(logit_if),
    contrast = colnames(contrast_matrix)[[i]],
    logit_difference = fit$coefficients[, i],
    usage_difference = raw_fit$coefficients[, i],
    p_value = fit$p.value[, i],
    stringsAsFactors = FALSE
  )
}))
pair_result <- merge(pair_result, contrast_metadata, by = "contrast", sort = FALSE)
nonfinite_pairwise_p_values <- sum(!is.finite(pair_result$p_value))
pair_result$p_value[!is.finite(pair_result$p_value)] <- 1
pair_result$isoform_q_global <- p.adjust(pair_result$p_value, method = "BH")
pair_result$isoform_q_by <- p.adjust(pair_result$p_value, method = "BY")

simes <- function(p) {
  p <- sort(p[is.finite(p)])
  if (!length(p)) return(1)
  min(1, min(p * length(p) / seq_along(p)))
}
pair_result <- merge(pair_result, mapping, by = "isoform_id", all.x = TRUE, sort = FALSE)
gene_p <- aggregate(p_value ~ gene_id + gene_name, pair_result, simes)
names(gene_p)[names(gene_p) == "p_value"] <- "gene_simes_p"
gene_p$gene_q <- p.adjust(gene_p$gene_simes_p, method = "BH")
bonferroni_min <- function(p) {
  p <- p[is.finite(p)]
  if (!length(p)) return(1)
  min(1, min(p) * length(p))
}
gene_conservative <- aggregate(
  p_value ~ gene_id + gene_name,
  pair_result,
  bonferroni_min
)
names(gene_conservative)[names(gene_conservative) == "p_value"] <-
  "gene_bonferroni_p"
gene_conservative$gene_q_by <- p.adjust(
  gene_conservative$gene_bonferroni_p,
  method = "BY"
)
gene_p <- merge(
  gene_p,
  gene_conservative,
  by = c("gene_id", "gene_name"),
  sort = FALSE
)
pair_result <- merge(pair_result, gene_p, by = c("gene_id", "gene_name"), all.x = TRUE, sort = FALSE)
pair_lookup <- split(pair_result, pair_result$contrast)

region_stage <- data.frame()
for (target in regions) {
  others <- setdiff(regions, target)
  for (stage in stage_levels) {
    comparison <- lapply(others, function(other) {
      pair <- sort(c(target, other))
      contrast_name <- paste(pair[[1]], "vs", pair[[2]], stage_codes[[stage]], sep = "_")
      value <- pair_lookup[[contrast_name]]
      value <- value[match(rownames(group_means), value$isoform_id), ]
      orientation <- if (target == pair[[1]]) 1 else -1
      data.frame(
        effect = orientation * value$usage_difference,
        q = value$isoform_q_global,
        q_by = value$isoform_q_by
      )
    })
    target_mean <- group_means[, paste(target, stage_codes[[stage]], sep = "_")]
    other_mean_1 <- group_means[, paste(others[[1]], stage_codes[[stage]], sep = "_")]
    other_mean_2 <- group_means[, paste(others[[2]], stage_codes[[stage]], sep = "_")]
    state <- ifelse(
      is.finite(comparison[[1]]$q) &
        is.finite(comparison[[2]]$q) &
        comparison[[1]]$q < 0.05 &
        comparison[[2]]$q < 0.05 &
        abs(comparison[[1]]$effect) >= 0.10 &
        abs(comparison[[2]]$effect) >= 0.10 &
        abs(other_mean_1 - other_mean_2) < 0.10 &
        sign(comparison[[1]]$effect) == sign(comparison[[2]]$effect),
      sign(comparison[[1]]$effect),
      0
    )
    region_stage <- rbind(
      region_stage,
      data.frame(
        isoform_id = rownames(group_means),
        region = target,
        stage = stage,
        stage_index = match(stage, stage_levels),
        target_mean_if = target_mean,
        other_mean_if = 0.5 * (other_mean_1 + other_mean_2),
        target_difference = target_mean - 0.5 * (other_mean_1 + other_mean_2),
        target_vs_other_1 = comparison[[1]]$effect,
        target_vs_other_2 = comparison[[2]]$effect,
        other_region_difference = other_mean_1 - other_mean_2,
        max_pair_q = pmax(comparison[[1]]$q, comparison[[2]]$q, na.rm = TRUE),
        max_pair_q_by = pmax(
          comparison[[1]]$q_by,
          comparison[[2]]$q_by,
          na.rm = TRUE
        ),
        state = state,
        stringsAsFactors = FALSE
      )
    )
  }
}
region_stage <- merge(region_stage, mapping, by = "isoform_id", all.x = TRUE, sort = FALSE)
region_stage <- merge(region_stage, gene_p, by = c("gene_id", "gene_name"), all.x = TRUE, sort = FALSE)

find_episodes <- function(
  stage_data,
  q_threshold = 0.05,
  effect_threshold = 0.10,
  flank_threshold = 0.10,
  q_column = "max_pair_q",
  gene_q_column = "gene_q"
) {
  value <- stage_data
  value$state <- ifelse(
    is.finite(value[[q_column]]) &
      value[[q_column]] < q_threshold &
      value[[gene_q_column]] < 0.05 &
      abs(value$target_vs_other_1) >= effect_threshold &
      abs(value$target_vs_other_2) >= effect_threshold &
      abs(value$other_region_difference) < effect_threshold &
      sign(value$target_vs_other_1) == sign(value$target_vs_other_2),
    sign(value$target_vs_other_1),
    0
  )
  groups <- split(
    value,
    interaction(value$region, value$isoform_id, drop = TRUE)
  )
  answer <- data.frame()
  for (group in groups) {
    isoform <- group$isoform_id[[1]]
    region <- group$region[[1]]
    group <- group[order(group$stage_index), ]
    runs <- rle(group$state)
    starts <- cumsum(c(1L, head(runs$lengths, -1L)))
    ends <- cumsum(runs$lengths)
    for (run_index in which(runs$values != 0)) {
      start <- starts[[run_index]]
      end <- ends[[run_index]]
      if (start <= 1L || end >= length(stage_levels)) next
      before <- start - 1L
      after <- end + 1L
      flanking_similarity <- max(
        abs(group$target_vs_other_1[c(before, after)]),
        abs(group$target_vs_other_2[c(before, after)]),
        na.rm = TRUE
      )
      if (!is.finite(flanking_similarity) ||
          flanking_similarity >= flank_threshold) {
        next
      }
      episode_rows <- start:end
      answer <- rbind(
        answer,
        data.frame(
          isoform_id = isoform,
          gene_id = group$gene_id[[1]],
          gene_name = group$gene_name[[1]],
          region = region,
          start_stage = group$stage[[start]],
          end_stage = group$stage[[end]],
          n_stages = end - start + 1L,
          direction = if (runs$values[[run_index]] > 0) "higher" else "lower",
          max_abs_usage_difference = max(abs(group$target_difference[episode_rows])),
          mean_abs_usage_difference = mean(abs(group$target_difference[episode_rows])),
          worst_pair_q = max(group[[q_column]][episode_rows], na.rm = TRUE),
          flanking_max_abs_difference = flanking_similarity,
          gene_q = group[[gene_q_column]][[1]],
          stringsAsFactors = FALSE
        )
      )
    }
  }
  answer[
    order(answer$worst_pair_q, -answer$max_abs_usage_difference, -answer$n_stages),
    ,
    drop = FALSE
  ]
}

episodes <- find_episodes(region_stage)

replicate_separation <- function(
  isoform_id,
  region,
  start_stage,
  end_stage,
  direction
) {
  episode_stages <- stage_levels[
    match(start_stage, stage_levels):match(end_stage, stage_levels)
  ]
  separations <- vapply(episode_stages, function(stage) {
    target_columns <- metadata$region == region &
      as.character(metadata$stage) == stage
    other_columns <- metadata$region != region &
      as.character(metadata$stage) == stage
    target_values <- isoform_fraction[isoform_id, target_columns]
    other_values <- isoform_fraction[isoform_id, other_columns]
    if (direction == "higher") {
      min(target_values) - max(other_values)
    } else {
      min(other_values) - max(target_values)
    }
  }, numeric(1))
  min(separations)
}
episodes$replicate_separation <- mapply(
  replicate_separation,
  episodes$isoform_id,
  episodes$region,
  episodes$start_stage,
  episodes$end_stage,
  episodes$direction
)
episodes$replicate_consistent <- episodes$replicate_separation > 0
statistical_episodes_before_replicate_check <- nrow(episodes)
episodes_with_incomplete_replicate_check <- sum(is.na(episodes$replicate_consistent))
episodes <- episodes[episodes$replicate_consistent %in% TRUE, , drop = FALSE]
episodes$episode_key <- paste(
  episodes$gene_id,
  episodes$region,
  episodes$start_stage,
  episodes$end_stage,
  sep = "__"
)
direction_count <- tapply(
  episodes$direction,
  episodes$episode_key,
  function(value) length(unique(value))
)
episodes$reciprocal_exchange <- direction_count[episodes$episode_key] > 1L
episodes$episode_key <- NULL

conservative_episodes <- find_episodes(
  region_stage,
  q_column = "max_pair_q_by",
  gene_q_column = "gene_q_by"
)
conservative_episodes$replicate_separation <- mapply(
  replicate_separation,
  conservative_episodes$isoform_id,
  conservative_episodes$region,
  conservative_episodes$start_stage,
  conservative_episodes$end_stage,
  conservative_episodes$direction
)
conservative_episodes$replicate_consistent <-
  conservative_episodes$replicate_separation > 0
conservative_episodes <- conservative_episodes[
  conservative_episodes$replicate_consistent %in% TRUE,
  ,
  drop = FALSE
]
conservative_key <- paste(
  conservative_episodes$gene_id,
  conservative_episodes$region,
  conservative_episodes$start_stage,
  conservative_episodes$end_stage,
  sep = "__"
)
conservative_direction_count <- tapply(
  conservative_episodes$direction,
  conservative_key,
  function(value) length(unique(value))
)
conservative_episodes$reciprocal_exchange <-
  conservative_direction_count[conservative_key] > 1L

episode_summary <- episodes |>
  dplyr::count(region, start_stage, end_stage, direction, name = "episodes") |>
  dplyr::arrange(match(start_stage, stage_levels), match(end_stage, stage_levels), direction)

sensitivity_settings <- data.frame(
  setting = c(
    "primary", "q_0.01", "effect_0.15", "effect_0.20",
    "flank_0.05", "flank_0.15"
  ),
  q_threshold = c(0.05, 0.01, 0.05, 0.05, 0.05, 0.05),
  effect_threshold = c(0.10, 0.10, 0.15, 0.20, 0.10, 0.10),
  flank_threshold = c(0.10, 0.10, 0.10, 0.10, 0.05, 0.15)
)
sensitivity <- bind_rows(lapply(seq_len(nrow(sensitivity_settings)), function(index) {
  setting <- sensitivity_settings[index, ]
  value <- find_episodes(
    region_stage,
    q_threshold = setting$q_threshold,
    effect_threshold = setting$effect_threshold,
    flank_threshold = setting$flank_threshold
  )
  if (nrow(value)) {
    value$replicate_separation <- mapply(
      replicate_separation,
      value$isoform_id,
      value$region,
      value$start_stage,
      value$end_stage,
      value$direction
    )
    value <- value[value$replicate_separation > 0 &
      !is.na(value$replicate_separation), , drop = FALSE]
  }
  data.frame(
    setting = setting$setting,
    q_threshold = setting$q_threshold,
    effect_threshold = setting$effect_threshold,
    flank_threshold = setting$flank_threshold,
    episodes = nrow(value),
    genes = length(unique(value$gene_id)),
    e15_5_single_stage_episodes = sum(
      value$start_stage == "15.5" & value$end_stage == "15.5"
    ),
    non_e15_5_episodes = sum(
      !(value$start_stage == "15.5" & value$end_stage == "15.5")
    )
  )
}))

reciprocal_groups <- episodes |>
  filter(reciprocal_exchange, replicate_consistent) |>
  group_by(gene_id, gene_name, region, start_stage, end_stage) |>
  summarise(
    isoforms = n(),
    max_worst_pair_q = max(worst_pair_q),
    max_abs_usage_difference = max(max_abs_usage_difference),
    min_replicate_separation = min(replicate_separation),
    .groups = "drop"
  ) |>
  arrange(max_worst_pair_q, desc(max_abs_usage_difference))
top_candidates <- reciprocal_groups |>
  distinct(gene_id, .keep_all = TRUE) |>
  slice_head(n = 6)

conservative_reciprocal_groups <- conservative_episodes |>
  filter(reciprocal_exchange, replicate_consistent) |>
  group_by(gene_id, gene_name, region, start_stage, end_stage) |>
  summarise(
    isoforms = n(),
    max_worst_pair_q = max(worst_pair_q),
    max_abs_usage_difference = max(max_abs_usage_difference),
    min_replicate_separation = min(replicate_separation),
    .groups = "drop"
  ) |>
  arrange(max_worst_pair_q, desc(max_abs_usage_difference))
conservative_top_candidates <- conservative_reciprocal_groups |>
  distinct(gene_id, .keep_all = TRUE) |>
  slice_head(n = 6)

episode_audit_row <- function(
  value,
  method,
  pairwise_adjustment,
  within_gene_aggregation,
  across_gene_adjustment,
  leading_genes
) {
  e15 <- value$start_stage == "15.5" & value$end_stage == "15.5"
  data.frame(
    method = method,
    pairwise_adjustment = pairwise_adjustment,
    within_gene_aggregation = within_gene_aggregation,
    across_gene_adjustment = across_gene_adjustment,
    replicate_separated_episodes = nrow(value),
    genes = length(unique(value$gene_id)),
    single_stage_e15_5_episodes = sum(e15),
    non_e15_5_episodes = sum(!e15),
    reciprocal_exchange_episodes = sum(value$reciprocal_exchange),
    reciprocal_exchange_genes = length(unique(
      value$gene_id[value$reciprocal_exchange]
    )),
    leading_six_genes = paste(leading_genes, collapse = ";"),
    stringsAsFactors = FALSE
  )
}
dependence_sensitivity <- bind_rows(
  episode_audit_row(
    episodes,
    "primary",
    "Benjamini-Hochberg over 300408 pairwise tests",
    "Simes over all pairwise tests per gene",
    "Benjamini-Hochberg over 4577 genes",
    top_candidates$gene_id
  ),
  episode_audit_row(
    conservative_episodes,
    "arbitrary-dependence sensitivity",
    "Benjamini-Yekutieli over 300408 pairwise tests",
    "Bonferroni minimum-p over all pairwise tests per gene",
    "Benjamini-Yekutieli over 4577 genes",
    conservative_top_candidates$gene_id
  )
)

episode_plot_data <- episodes |>
  dplyr::count(start_stage, name = "episodes") |>
  mutate(
    start_stage = factor(start_stage, levels = stage_levels),
    e15_5 = start_stage == "15.5"
  )
p_episode_counts <- ggplot(
  episode_plot_data,
  aes(start_stage, episodes, fill = e15_5)
) +
  geom_col(width = 0.72) +
  geom_text(
    aes(start_stage, episodes, label = episodes),
    inherit.aes = FALSE,
    vjust = -0.35,
    size = 3.0
  ) +
  scale_fill_manual(values = c(`FALSE` = "#6BAED6", `TRUE` = "#D55E00")) +
  scale_x_discrete(labels = stage_labels) +
  scale_y_sqrt(expand = expansion(mult = c(0, 0.10))) +
  labs(
    title = "A  Statistically calibrated transient regional episodes",
    subtitle = "Episode start stage; square-root y scale makes non-E15.5 calls visible",
    x = "Developmental stage", y = "Isoform episodes (square-root scale)"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.title = element_text(face = "bold")
  )

candidate_episodes <- episodes |>
  semi_join(
    top_candidates,
    by = c("gene_id", "gene_name", "region", "start_stage", "end_stage")
  )
candidate_isoforms <- candidate_episodes |>
  dplyr::select(
    gene_id, gene_name, region, start_stage, end_stage, isoform_id,
    direction, worst_pair_q, max_abs_usage_difference, replicate_separation
  ) |>
  arrange(gene_name, direction, worst_pair_q)
candidate_ids <- unique(candidate_episodes$isoform_id)
candidate_if <- as.data.frame(
  isoform_fraction[candidate_ids, , drop = FALSE],
  check.names = FALSE
) |>
  mutate(isoform_id = rownames(isoform_fraction[candidate_ids, , drop = FALSE])) |>
  pivot_longer(-isoform_id, names_to = "sample_id", values_to = "isoform_fraction") |>
  left_join(
    metadata[, c("sample_id", "region", "stage", "replicate")],
    by = "sample_id"
  ) |>
  left_join(mapping[, c("isoform_id", "gene_name")], by = "isoform_id") |>
  mutate(stage = factor(as.character(stage), levels = stage_levels))
candidate_means <- candidate_if |>
  left_join(
    candidate_isoforms |> distinct(gene_name, isoform_id, direction),
    by = c("gene_name", "isoform_id")
  ) |>
  mutate(region_group = ifelse(region == "Midbrain", "Midbrain", "Forebrain/hindbrain mean")) |>
  group_by(gene_name, isoform_id, direction, region_group, stage) |>
  summarise(isoform_fraction = mean(isoform_fraction), .groups = "drop")

p_candidates <- ggplot(
  candidate_means,
  aes(
    stage,
    isoform_fraction,
    colour = direction,
    linetype = region_group,
    group = interaction(isoform_id, region_group)
  )
) +
  annotate("rect", xmin = 5.5, xmax = 6.5, ymin = -Inf, ymax = Inf,
           fill = "#F0E442", alpha = 0.12) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.35) +
  facet_wrap(~ gene_name, scales = "free_y", ncol = 3) +
  scale_colour_manual(
    values = c(higher = "#D55E00", lower = "#0072B2"),
    labels = c(higher = "Higher at episode", lower = "Lower at episode")
  ) +
  scale_x_discrete(labels = stage_labels) +
  scale_linetype_manual(values = c("Midbrain" = "solid", "Forebrain/hindbrain mean" = "22")) +
  labs(
    title = "B  Top replicate-consistent reciprocal exchanges",
    subtitle = "Higher/lower denote the midbrain fraction at E15.5; RefSeq accessions are reported with the figure",
    x = "Developmental stage", y = "Isoform fraction",
    colour = NULL, linetype = NULL
  ) +
  theme_minimal(base_size = 9.5) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

figure <- p_episode_counts / p_candidates + plot_layout(heights = c(0.65, 1.4))
ggsave(
  file.path(figure_dir, "figure6_transient_regional_episodes.pdf"),
  figure,
  width = 7.2,
  height = 8.2,
  units = "in",
  device = cairo_pdf
)
ggsave(
  file.path(figure_dir, "figure6_transient_regional_episodes.png"),
  figure,
  width = 7.2,
  height = 8.2,
  units = "in",
  dpi = 300,
  bg = "white"
)

significant_pair_tests <- pair_result |>
  filter(
    is.finite(isoform_q_global),
    isoform_q_global < 0.05,
    abs(usage_difference) >= 0.10
  ) |>
  arrange(isoform_q_global, desc(abs(usage_difference)))
significant_stage_calls <- region_stage |>
  filter(state != 0) |>
  arrange(max_pair_q, desc(abs(target_difference)))

pair_result_release <- pair_result |>
  dplyr::select(
    isoform_id, gene_id, gene_name, contrast, region_1, region_2, stage,
    logit_difference, usage_difference, p_value, isoform_q_global,
    isoform_q_by, gene_simes_p, gene_q, gene_bonferroni_p, gene_q_by
  ) |>
  arrange(isoform_id, region_1, region_2, match(stage, stage_levels))

region_stage_release <- region_stage |>
  dplyr::select(
    isoform_id, gene_id, gene_name, region, stage, stage_index,
    target_mean_if, other_mean_if, target_difference,
    target_vs_other_1, target_vs_other_2, other_region_difference,
    max_pair_q, max_pair_q_by, gene_simes_p, gene_q,
    gene_bonferroni_p, gene_q_by
  ) |>
  arrange(isoform_id, region, stage_index)

isoform_fraction_release <- data.frame(
  isoform_id = rownames(isoform_fraction),
  isoform_fraction,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

diagnostics <- data.frame(
  metric = c(
    "common_isoforms_before_filtering",
    "isoforms_after_expression_and_multi_isoform_filtering",
    "multi_isoform_genes_tested",
    "samples",
    "region_stage_cells",
    "pairwise_tests",
    "nonfinite_pairwise_p_values_replaced_by_bh_as_one",
    "statistical_episodes_before_replicate_check",
    "episodes_with_incomplete_replicate_check_excluded",
    "primary_episodes",
    "primary_episode_genes",
    "replicate_consistent_episodes",
    "reciprocal_exchange_episodes",
    "replicate_consistent_reciprocal_genes",
    "single_stage_e15_5_episodes",
    "non_e15_5_episodes"
  ),
  value = c(
    length(common_isoforms),
    nrow(mapping),
    length(unique(mapping$gene_id)),
    ncol(isoform_fraction),
    nlevels(metadata$group),
    nrow(pair_result),
    nonfinite_pairwise_p_values,
    statistical_episodes_before_replicate_check,
    episodes_with_incomplete_replicate_check,
    nrow(episodes),
    length(unique(episodes$gene_id)),
    sum(episodes$replicate_consistent),
    sum(episodes$reciprocal_exchange),
    length(unique(
      episodes$gene_id[episodes$replicate_consistent & episodes$reciprocal_exchange]
    )),
    sum(episodes$start_stage == "15.5" & episodes$end_stage == "15.5"),
    sum(!(episodes$start_stage == "15.5" & episodes$end_stage == "15.5"))
  )
)

write.csv(
  significant_pair_tests,
  file.path(data_dir, "transient_regional_pair_tests_significant.csv"),
  row.names = FALSE
)
write.csv(
  pair_result_release,
  file.path(data_dir, "transient_regional_pair_tests_all.csv"),
  row.names = FALSE
)
write.csv(
  significant_stage_calls,
  file.path(data_dir, "transient_regional_stage_calls.csv"),
  row.names = FALSE
)
write.csv(
  region_stage_release,
  file.path(data_dir, "transient_regional_stage_evaluations_all.csv"),
  row.names = FALSE
)
write.csv(
  isoform_fraction_release,
  file.path(data_dir, "transient_regional_filtered_isoform_fractions.csv"),
  row.names = FALSE
)
write.csv(
  episodes,
  file.path(data_dir, "transient_regional_isoform_episodes.csv"),
  row.names = FALSE
)
write.csv(
  conservative_episodes,
  file.path(data_dir, "transient_regional_isoform_episodes_conservative.csv"),
  row.names = FALSE
)
write.csv(
  episode_summary,
  file.path(table_dir, "transient_regional_episode_summary.csv"),
  row.names = FALSE
)
write.csv(
  sensitivity,
  file.path(table_dir, "transient_regional_sensitivity.csv"),
  row.names = FALSE
)
write.csv(
  dependence_sensitivity,
  file.path(table_dir, "transient_regional_dependence_sensitivity.csv"),
  row.names = FALSE
)
write.csv(
  top_candidates,
  file.path(table_dir, "transient_regional_top_candidates.csv"),
  row.names = FALSE
)
write.csv(
  candidate_isoforms,
  file.path(table_dir, "transient_regional_top_candidate_isoforms.csv"),
  row.names = FALSE
)
write.csv(
  candidate_means,
  file.path(table_dir, "transient_regional_top_candidate_trajectories.csv"),
  row.names = FALSE
)
write.csv(
  diagnostics,
  file.path(table_dir, "transient_regional_scan_diagnostics.csv"),
  row.names = FALSE
)

cat("Transient regional scan complete.\n")
print(diagnostics)
cat("\nThreshold sensitivity:\n")
print(sensitivity)
cat("\nDependence sensitivity:\n")
print(dependence_sensitivity)
cat("\nTop reciprocal candidates:\n")
print(top_candidates)
