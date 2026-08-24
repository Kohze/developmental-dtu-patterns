#!/usr/bin/env Rscript

# Reconstruct replicate-level methylation summaries and transcript usage from
# the archived thesis inputs. This script does not quantile-normalize or pair
# independent WGBS and RNA replicates.

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(data.table)
})

arguments <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", grep("^--file=", arguments, value = TRUE)[1])
paper_dir <- normalizePath(dirname(script_file), winslash = "/", mustWork = TRUE)
root <- normalizePath(file.path(paper_dir, "..", ".."), winslash = "/", mustWork = TRUE)
cli_arguments <- commandArgs(trailingOnly = TRUE)
input_dir_argument <- if (length(cli_arguments)) cli_arguments[[1]] else ""
input_dir_environment <- Sys.getenv("METHYLATION_INPUT_DIR", unset = "")
input_dir <- if (nzchar(input_dir_argument)) {
  normalizePath(input_dir_argument, winslash = "/", mustWork = TRUE)
} else if (nzchar(input_dir_environment)) {
  normalizePath(input_dir_environment, winslash = "/", mustWork = TRUE)
} else {
  root
}
results_dir <- file.path(paper_dir, "results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

sample_manifest <- rbindlist(list(
  data.table(
    tissue = "ForeBrain",
    stage = rep(c("10.5", "11.5", "12.5", "13.5", "14.5", "15.5", "16.5", "0"), each = 2),
    replicate = rep(1:2, 8),
    sample_id = c(
      "ENCFF288JMC", "ENCFF461PTM", "ENCFF890XTY", "ENCFF232QBX",
      "ENCFF535JWX", "ENCFF667LWJ", "ENCFF673PJA", "ENCFF320IBB",
      "ENCFF667IZN", "ENCFF948AII", "ENCFF331SNN", "ENCFF566THP",
      "ENCFF457SPJ", "ENCFF480MYA", "ENCFF486VKI", "ENCFF784IJC"
    )
  ),
  data.table(
    tissue = "MidBrain",
    stage = rep(c("10.5", "11.5", "12.5", "13.5", "14.5", "15.5", "16.5", "0"), each = 2),
    replicate = rep(1:2, 8),
    sample_id = c(
      "ENCFF192LFV", "ENCFF724ZCI", "ENCFF943IBD", "ENCFF580ZKC",
      "ENCFF910CHR", "ENCFF147HUB", "ENCFF376OQK", "ENCFF323DLV",
      "ENCFF496WRJ", "ENCFF499EUN", "ENCFF354OEF", "ENCFF678MHU",
      "ENCFF451ZEE", "ENCFF529OUU", "ENCFF415LLK", "ENCFF733TCH"
    )
  )
))
fwrite(sample_manifest, file.path(results_dir, "wgbs_sample_manifest.csv"))

environment <- new.env(parent = emptyenv())
load(file.path(input_dir, "methylation_subset_grl.RData"), envir = environment)
methylation_samples <- environment$methylation_subset_grl
stopifnot(all(sample_manifest$sample_id %in% names(methylation_samples)))

summarize_one_sample <- function(regions, sample_ranges, metadata) {
  # CpG methylation is not directional with respect to gene annotation. The
  # original findOverlaps call inherited strand sensitivity; the primary
  # reconstruction explicitly ignores strand and records the original
  # strand-sensitive estimate as a sensitivity column.
  summarize_overlap <- function(ignore_strand) {
    overlap <- findOverlaps(regions, sample_ranges, ignore.strand = ignore_strand)
    if (!length(overlap)) {
      return(data.table(region_index = integer(), n_cpg = integer(),
                        mean_unweighted = numeric(), mean_weighted = numeric(),
                        total_coverage = numeric()))
    }
    query <- queryHits(overlap)
    subject <- subjectHits(overlap)
    values <- as.numeric(mcols(sample_ranges)$methylation_level[subject])
    coverage <- as.numeric(mcols(sample_ranges)$Coverage[subject])
    data.table(query = query, values = values, coverage = coverage)[
      is.finite(values) & is.finite(coverage) & coverage > 0,
      .(
        n_cpg = .N,
        mean_unweighted = mean(values),
        mean_weighted = weighted.mean(values, coverage),
        total_coverage = sum(coverage)
      ),
      by = .(region_index = query)
    ]
  }

  primary <- summarize_overlap(TRUE)
  original <- summarize_overlap(FALSE)
  setnames(
    original,
    c("n_cpg", "mean_unweighted", "mean_weighted", "total_coverage"),
    paste0(c("n_cpg", "mean_unweighted", "mean_weighted", "total_coverage"),
           "_strand_sensitive")
  )
  output <- data.table(region_index = seq_along(regions))
  output <- merge(output, primary, by = "region_index", all.x = TRUE)
  output <- merge(output, original, by = "region_index", all.x = TRUE)
  output[, `:=`(
    tissue = metadata$tissue,
    stage = metadata$stage,
    replicate = metadata$replicate,
    sample_id = metadata$sample_id
  )]
  output
}

extract_tissue_methylation <- function(tissue_name, rds_file) {
  message("Reconstructing methylation for ", tissue_name)
  archived <- readRDS(file.path(input_dir, rds_file))
  region_metadata <- data.table(
    region_index = seq_along(archived),
    region_key = names(archived),
    gene_id = vapply(archived, function(item) item$isoform_id, character(1)),
    region = vapply(archived, function(item) item$region, character(1)),
    original_width = vapply(archived, function(item) width(item$range), integer(1)),
    analysed_width = vapply(archived, function(item) width(item$extended_range), integer(1))
  )
  regions <- GRangesList(lapply(archived, function(item) item$extended_range))
  regions <- unlist(regions, use.names = FALSE)
  stopifnot(length(regions) == nrow(region_metadata))

  manifest <- sample_manifest[sample_manifest$tissue == tissue_name]
  sample_results <- lapply(seq_len(nrow(manifest)), function(index) {
    metadata <- manifest[index]
    message("  ", metadata$sample_id, " (", metadata$stage, "_", metadata$replicate, ")")
    summarize_one_sample(
      regions,
      methylation_samples[[metadata$sample_id]],
      metadata
    )
  })
  output <- rbindlist(sample_results)
  merge(output, region_metadata, by = "region_index", all.x = TRUE)
}

forebrain_methylation <- extract_tissue_methylation(
  "ForeBrain", "forebrain_isoforms_data_with_dmr.rds"
)
midbrain_methylation <- extract_tissue_methylation(
  "MidBrain", "midbrain_isoforms_data_with_dmr.rds"
)
replicate_methylation <- rbindlist(list(forebrain_methylation, midbrain_methylation))
setcolorder(
  replicate_methylation,
  c("tissue", "gene_id", "region", "region_key", "stage", "replicate",
    "sample_id", "original_width", "analysed_width", "n_cpg",
    "mean_unweighted", "mean_weighted", "total_coverage",
    "n_cpg_strand_sensitive", "mean_unweighted_strand_sensitive",
    "mean_weighted_strand_sensitive", "total_coverage_strand_sensitive",
    "region_index")
)
fwrite(
  replicate_methylation,
  file.path(results_dir, "replicate_level_methylation.csv")
)

extract_expression <- function(tissue, rdata_file, object_name) {
  message("Extracting transcript usage for ", tissue)
  environment <- new.env()
  load(file.path(input_dir, rdata_file), envir = environment)
  object <- environment[[object_name]]

  usage <- as.data.table(object$isoformRepIF)
  expression <- as.data.table(object$isoformRepExpression)
  long_usage <- melt(
    usage,
    id.vars = "isoform_id",
    variable.name = "rna_sample",
    value.name = "rep_if"
  )
  long_expression <- melt(
    expression,
    id.vars = "isoform_id",
    variable.name = "rna_sample",
    value.name = "expression"
  )
  output <- merge(
    long_usage, long_expression,
    by = c("isoform_id", "rna_sample"),
    all = TRUE
  )
  output[, stage := sub("_.*$", "", rna_sample)]
  output[, replicate := as.integer(sub("^.*_", "", rna_sample))]
  output[, tissue := tissue]
  output[]
}

replicate_expression <- rbindlist(list(
  extract_expression("ForeBrain", "foreBrain.RData", "combinedForebrain"),
  extract_expression("MidBrain", "isoform_final.RData", "combinedMitbrain")
))
setcolorder(
  replicate_expression,
  c("tissue", "isoform_id", "stage", "replicate", "rna_sample",
    "rep_if", "expression")
)
fwrite(
  replicate_expression,
  file.path(results_dir, "replicate_level_isoform_usage.csv")
)

session <- capture.output(sessionInfo())
writeLines(session, file.path(results_dir, "reconstruction_session_info.txt"))
message("Reconstruction complete.")
