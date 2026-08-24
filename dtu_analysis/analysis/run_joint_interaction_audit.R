#!/usr/bin/env Rscript

# Focused joint region-by-stage DTU audit for the E15.5-centred midbrain claim.
# Usage:
#   Rscript analysis/run_joint_interaction_audit.R [directory containing the two .RData inputs]

options(stringsAsFactors = FALSE, width = 140)
set.seed(20260730)

suppressPackageStartupMessages({
  library(satuRn)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(edgeR)
  library(DEXSeq)
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
  Forebrain = switchListForebrain,
  Hindbrain = switchListHindbrain,
  Midbrain = switchListMidbrain
)

extract_counts <- function(x) {
  value <- x$isoformCountMatrix
  stopifnot(!anyDuplicated(value$isoform_id))
  matrix_value <- as.matrix(value[, -1, drop = FALSE])
  storage.mode(matrix_value) <- "numeric"
  rownames(matrix_value) <- value$isoform_id
  matrix_value
}

extract_mapping <- function(x) {
  mapping <- unique(x$isoformFeatures[, c("isoform_id", "gene_id", "gene_name")])
  mapping_count <- aggregate(
    gene_id ~ isoform_id,
    data = mapping,
    FUN = function(value) length(unique(value))
  )
  if (any(mapping_count$gene_id != 1L)) {
    stop("An isoform maps to more than one gene in an archived object.")
  }
  mapping[!duplicated(mapping$isoform_id), , drop = FALSE]
}

count_matrices <- lapply(objects, extract_counts)
mapping_tables <- lapply(objects, extract_mapping)
common_isoforms <- Reduce(intersect, lapply(count_matrices, rownames))

reference_mapping <- mapping_tables$Forebrain[
  match(common_isoforms, mapping_tables$Forebrain$isoform_id),
  c("isoform_id", "gene_id", "gene_name")
]
if (anyNA(reference_mapping$gene_id)) {
  stop("The forebrain mapping is incomplete for common isoforms.")
}
for (region in c("Hindbrain", "Midbrain")) {
  comparison <- mapping_tables[[region]][
    match(common_isoforms, mapping_tables[[region]]$isoform_id),
    ,
    drop = FALSE
  ]
  if (!identical(
    as.character(reference_mapping$gene_id),
    as.character(comparison$gene_id)
  )) {
    stop("Cross-region isoform-to-gene mappings do not agree for ", region, ".")
  }
}

combined_counts <- do.call(cbind, lapply(names(count_matrices), function(region) {
  value <- count_matrices[[region]][common_isoforms, , drop = FALSE]
  colnames(value) <- paste(region, colnames(value), sep = "__")
  value
}))

sample_parts <- strsplit(colnames(combined_counts), "__", fixed = TRUE)
sample_metadata <- data.frame(
  sample_id = colnames(combined_counts),
  region = vapply(sample_parts, `[[`, character(1), 1),
  stage_replicate = vapply(sample_parts, `[[`, character(1), 2),
  stringsAsFactors = FALSE
)
sample_metadata$stage <- sub("_[12]$", "", sample_metadata$stage_replicate)
sample_metadata$replicate <- sub("^.*_", "", sample_metadata$stage_replicate)
stage_levels <- c("10.5", "11.5", "12.5", "13.5", "14.5", "15.5", "16.5", "0")
stage_codes <- c(
  "10.5" = "E10_5", "11.5" = "E11_5", "12.5" = "E12_5",
  "13.5" = "E13_5", "14.5" = "E14_5", "15.5" = "E15_5",
  "16.5" = "E16_5", "0" = "P0"
)
sample_metadata$stage <- factor(sample_metadata$stage, levels = stage_levels)
sample_metadata$region <- factor(
  sample_metadata$region,
  levels = c("Forebrain", "Hindbrain", "Midbrain")
)
sample_metadata$stage_code <- unname(stage_codes[as.character(sample_metadata$stage)])
group_levels <- unlist(lapply(
  levels(sample_metadata$region),
  function(region) paste(region, unname(stage_codes[stage_levels]), sep = "_")
))
sample_metadata$group <- factor(
  paste(sample_metadata$region, sample_metadata$stage_code, sep = "_"),
  levels = group_levels
)
rownames(sample_metadata) <- sample_metadata$sample_id

design <- model.matrix(~ 0 + group, sample_metadata)
colnames(design) <- levels(sample_metadata$group)
if (qr(design)$rank != ncol(design)) {
  stop("The joint region-by-stage design is not full rank.")
}
cell_counts <- table(sample_metadata$group)
if (any(cell_counts != 2L)) {
  stop("Expected exactly two samples in every region-stage cell.")
}

expression_filter <- filterByExpr(
  combined_counts,
  group = sample_metadata$group,
  min.count = 10,
  min.total.count = 30,
  large.n = 20,
  min.prop = 0.7
)
filtered_mapping <- reference_mapping[expression_filter, , drop = FALSE]
multi_isoform_gene <- duplicated(filtered_mapping$gene_id) |
  duplicated(filtered_mapping$gene_id, fromLast = TRUE)
analysis_mapping <- filtered_mapping[multi_isoform_gene, , drop = FALSE]
analysis_counts <- combined_counts[
  analysis_mapping$isoform_id,
  ,
  drop = FALSE
]
rownames(analysis_mapping) <- analysis_mapping$isoform_id

if (!identical(rownames(analysis_counts), rownames(analysis_mapping))) {
  stop("Count and mapping rows are not aligned.")
}

sum_exp <- SummarizedExperiment(
  assays = list(counts = analysis_counts),
  colData = DataFrame(sample_metadata),
  rowData = DataFrame(analysis_mapping)
)
metadata(sum_exp)$formula <- ~ 0 + group

fit_warning_messages <- character()
sum_exp <- withCallingHandlers(
  fitDTU(
    object = sum_exp,
    formula = ~ 0 + group,
    parallel = FALSE,
    verbose = TRUE
  ),
  warning = function(warning_condition) {
    fit_warning_messages <<- c(
      fit_warning_messages,
      conditionMessage(warning_condition)
    )
    invokeRestart("muffleWarning")
  }
)
fit_model_types <- vapply(
  rowData(sum_exp)[["fitDTUModels"]],
  function(model) model@type,
  character(1)
)

contrast_name <- "E15_5_midbrain_curvature"
contrast <- matrix(
  0,
  nrow = ncol(design),
  ncol = 1,
  dimnames = list(colnames(design), contrast_name)
)
contrast[c("Midbrain_E14_5", "Midbrain_E15_5", "Midbrain_E16_5"), 1] <-
  c(-0.5, 1, -0.5)
contrast[c("Forebrain_E14_5", "Forebrain_E15_5", "Forebrain_E16_5"), 1] <-
  c(0.25, -0.5, 0.25)
contrast[c("Hindbrain_E14_5", "Hindbrain_E15_5", "Hindbrain_E16_5"), 1] <-
  c(0.25, -0.5, 0.25)
if (abs(sum(contrast)) > 1e-12) {
  stop("The interaction contrast does not sum to zero.")
}

test_warning_messages <- character()
# satuRn 1.12.0 draws its empirical-null diagnostics even when both diagplot
# arguments are FALSE. Route those plots to a null device so a reconstruction
# does not leave an undocumented Rplots.pdf artifact.
grDevices::pdf(file = NULL)
sum_exp <- tryCatch(
  withCallingHandlers(
    testDTU(
      object = sum_exp,
      contrasts = contrast,
      diagplot1 = FALSE,
      diagplot2 = FALSE,
      sort = FALSE
    ),
    warning = function(warning_condition) {
      test_warning_messages <<- c(
        test_warning_messages,
        conditionMessage(warning_condition)
      )
      invokeRestart("muffleWarning")
    }
  ),
  finally = grDevices::dev.off()
)

result_field <- paste0("fitDTUResult_", contrast_name)
isoform_result <- as.data.frame(rowData(sum_exp)[[result_field]])
isoform_result <- cbind(
  as.data.frame(rowData(sum_exp)[, c("isoform_id", "gene_id", "gene_name")]),
  isoform_result
)
invalid_empirical_p_values <- sum(!is.finite(isoform_result$empirical_pval))
invalid_empirical_fdr_values <- sum(!is.finite(isoform_result$empirical_FDR))
isoform_result$inference_valid <- is.finite(isoform_result$empirical_pval) &
  is.finite(isoform_result$empirical_FDR)
isoform_result$empirical_pval[!is.finite(isoform_result$empirical_pval)] <- 1
isoform_result$empirical_FDR[!is.finite(isoform_result$empirical_FDR)] <- 1

gene_factor <- factor(isoform_result$gene_id)
gene_split <- split(seq_along(gene_factor), gene_factor)
gene_min_p <- vapply(
  gene_split,
  function(index) min(isoform_result$empirical_pval[index]),
  numeric(1)
)
gene_min_p[!is.finite(gene_min_p)] <- 1
theta <- unique(sort(gene_min_p))
gene_q_lookup <- DEXSeq:::perGeneQValueExact(gene_min_p, theta, gene_split)
gene_q <- pmin(1, gene_q_lookup[match(gene_min_p, theta)])
names(gene_q) <- names(gene_split)

best_index <- vapply(
  gene_split,
  function(index) index[which.min(isoform_result$empirical_pval[index])],
  integer(1)
)
gene_result <- data.frame(
  gene_id = names(gene_split),
  gene_name = isoform_result$gene_name[best_index],
  n_tested_isoforms = lengths(gene_split),
  gene_min_empirical_p = unname(gene_min_p),
  gene_q_value = unname(gene_q),
  best_isoform = isoform_result$isoform_id[best_index],
  best_isoform_estimate = isoform_result$estimates[best_index],
  best_isoform_empirical_p = isoform_result$empirical_pval[best_index],
  best_isoform_empirical_FDR = isoform_result$empirical_FDR[best_index],
  stringsAsFactors = FALSE
)

membership_path <- file.path(data_dir, "dtu_gene_membership.csv")
if (file.exists(membership_path)) {
  membership <- read.csv(membership_path, check.names = FALSE)
  gene_result <- merge(
    gene_result,
    membership[, c("gene_id", "Forebrain", "Hindbrain", "Midbrain", "membership")],
    by = "gene_id",
    all.x = TRUE,
    sort = FALSE
  )
}
gene_result <- gene_result[
  order(gene_result$gene_q_value, gene_result$gene_min_empirical_p),
  ,
  drop = FALSE
]

contrast_table <- data.frame(
  group = rownames(contrast),
  coefficient = as.numeric(contrast[, 1]),
  samples = as.integer(cell_counts[rownames(contrast)]),
  stringsAsFactors = FALSE
)

summarise_warnings <- function(messages, phase) {
  if (!length(messages)) {
    return(data.frame(
      phase = character(),
      message = character(),
      count = integer(),
      stringsAsFactors = FALSE
    ))
  }
  warning_counts <- sort(table(messages), decreasing = TRUE)
  data.frame(
    phase = phase,
    message = names(warning_counts),
    count = as.integer(warning_counts),
    stringsAsFactors = FALSE
  )
}
warning_table <- rbind(
  summarise_warnings(fit_warning_messages, "fitDTU"),
  summarise_warnings(test_warning_messages, "testDTU")
)

significant_genes <- gene_result$gene_q_value < 0.05
genes_with_any_invalid_test <- unique(
  isoform_result$gene_id[!isoform_result$inference_valid]
)
genes_with_all_tests_invalid <- names(gene_split)[vapply(
  gene_split,
  function(index) all(!isoform_result$inference_valid[index]),
  logical(1)
)]
gene_result$any_nonfinite_transcript_test <-
  gene_result$gene_id %in% genes_with_any_invalid_test
gene_result$all_transcript_tests_nonfinite <-
  gene_result$gene_id %in% genes_with_all_tests_invalid
candidate_rows <- gene_result[gene_result$gene_id %in% c("Ppp2r3a", "Rtn2"), ]

if (any(significant_genes & gene_result$all_transcript_tests_nonfinite)) {
  stop("A gene with no finite transcript test was called significant.")
}
if (any(
  !isoform_result$inference_valid &
    isoform_result$gene_id %in% c("Ppp2r3a", "Rtn2")
)) {
  stop("A reported candidate has a nonfinite transcript test.")
}
primary_shared <- !is.na(gene_result$membership) &
  gene_result$membership == "Forebrain + Hindbrain + Midbrain"
primary_midbrain_only <- !is.na(gene_result$membership) &
  gene_result$membership == "Midbrain only"

metric <- c(
  "common_isoforms_before_filtering",
  "isoforms_passing_expression_filter",
  "isoforms_tested_in_multi_isoform_genes",
  "multi_isoform_genes_tested",
  "samples",
  "region_stage_cells",
  "design_rank",
  "design_residual_df",
  "noninteger_count_fraction",
  "fit_warning_count",
  "fit_nonconvergence_warning_count",
  "test_warning_count",
  "fit_error_models",
  "invalid_empirical_p_values_replaced_with_1",
  "invalid_empirical_fdr_values_replaced_with_1",
  "genes_with_any_nonfinite_transcript_test",
  "genes_with_all_transcript_tests_nonfinite",
  "significant_genes_with_any_nonfinite_transcript_test",
  "candidate_nonfinite_transcript_tests",
  "transcripts_empirical_fdr_lt_0_05",
  "genes_dexseq_style_per_gene_q_lt_0_05",
  "significant_genes_in_primary_shared_core",
  "significant_genes_in_primary_midbrain_only",
  "Ppp2r3a_gene_q_value",
  "Rtn2_gene_q_value"
)
value <- c(
  length(common_isoforms),
  sum(expression_filter),
  nrow(analysis_counts),
  length(unique(analysis_mapping$gene_id)),
  ncol(analysis_counts),
  nlevels(sample_metadata$group),
  qr(design)$rank,
  nrow(design) - qr(design)$rank,
  mean(abs(analysis_counts - round(analysis_counts)) > 1e-8),
  length(fit_warning_messages),
  sum(grepl("did not converge", fit_warning_messages, fixed = TRUE)),
  length(test_warning_messages),
  sum(fit_model_types == "fitError"),
  invalid_empirical_p_values,
  invalid_empirical_fdr_values,
  length(genes_with_any_invalid_test),
  length(genes_with_all_tests_invalid),
  sum(significant_genes & gene_result$gene_id %in% genes_with_any_invalid_test),
  sum(
    !isoform_result$inference_valid &
      isoform_result$gene_id %in% c("Ppp2r3a", "Rtn2")
  ),
  sum(isoform_result$empirical_FDR < 0.05, na.rm = TRUE),
  sum(significant_genes, na.rm = TRUE),
  sum(significant_genes & primary_shared, na.rm = TRUE),
  sum(significant_genes & primary_midbrain_only, na.rm = TRUE),
  candidate_rows$gene_q_value[match("Ppp2r3a", candidate_rows$gene_id)],
  candidate_rows$gene_q_value[match("Rtn2", candidate_rows$gene_id)]
)
summary_table <- data.frame(metric = metric, value = value)

write.csv(
  isoform_result,
  file.path(data_dir, "joint_interaction_isoforms.csv"),
  row.names = FALSE
)
write.csv(
  gene_result,
  file.path(data_dir, "joint_interaction_genes.csv"),
  row.names = FALSE
)
write.csv(
  summary_table,
  file.path(table_dir, "joint_interaction_summary.csv"),
  row.names = FALSE
)
write.csv(
  contrast_table,
  file.path(table_dir, "joint_interaction_contrast.csv"),
  row.names = FALSE
)
write.csv(
  warning_table,
  file.path(table_dir, "joint_interaction_warnings.csv"),
  row.names = FALSE
)

cat("Joint interaction audit complete.\n")
print(summary_table)
cat("\nTop joint-interaction genes:\n")
print(head(gene_result, 15))
cat("\nCandidate rows:\n")
print(candidate_rows)
