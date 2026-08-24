#!/usr/bin/env Rscript

# Audit the recoverable quantification and annotation structure in the archived
# workspaces without altering either input.
# Usage:
#   Rscript analysis/audit_archived_inputs.R [directory containing the two .RData inputs]

options(stringsAsFactors = FALSE, width = 120)

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

load_isolated <- function(filename) {
  environment <- new.env(parent = emptyenv())
  load(file.path(input_dir, filename), envir = environment)
  environment
}

forebrain_workspace <- load_isolated("foreBrain.RData")
midbrain_workspace <- load_isolated("isoform_final.RData")

refseq_profile <- function(ids) {
  base_ids <- sub("^((NM|NR)_[0-9]+)_[0-9]+$", "\\1", ids)
  list(
    all_match = all(grepl("^(NM|NR)_[0-9]+(_[0-9]+)?$", ids)),
    base_accessions = length(unique(base_ids)),
    suffix_instances = sum(ids != base_ids)
  )
}

import_options <- function(quantification) {
  paste(
    names(quantification$importOptions),
    unlist(quantification$importOptions, use.names = FALSE),
    sep = "=",
    collapse = "; "
  )
}

count_lineage <- function(quantification, switch_object) {
  imported <- quantification$counts
  retained <- switch_object$isoformCountMatrix
  retained_ids <- retained$isoform_id
  imported_index <- match(retained_ids, imported$isoform_id)
  sample_columns <- intersect(
    setdiff(colnames(retained), "isoform_id"),
    setdiff(colnames(imported), "isoform_id")
  )
  if (anyNA(imported_index) || !length(sample_columns)) {
    return(FALSE)
  }
  isTRUE(all.equal(
    as.matrix(retained[, sample_columns, drop = FALSE]),
    as.matrix(imported[imported_index, sample_columns, drop = FALSE]),
    check.attributes = FALSE,
    tolerance = 0
  ))
}

workspace_quantification <- list(
  Forebrain = forebrain_workspace$salmonQuant,
  Midbrain = midbrain_workspace$salmonQuant
)
shared_quantification_ids <- identical(
  workspace_quantification$Forebrain$abundance$isoform_id,
  workspace_quantification$Midbrain$abundance$isoform_id
)
shared_length_matrix <- identical(
  workspace_quantification$Forebrain$length,
  workspace_quantification$Midbrain$length
)
cross_workspace_evidence <- paste0(
  "Forebrain and midbrain salmonQuant objects have identical transcript-ID vectors=",
  shared_quantification_ids,
  " and identical length matrices=",
  shared_length_matrix,
  "; abundance and count matrices differ"
)

regions <- list(
  Forebrain = list(
    workspace = "foreBrain.RData",
    switch_object = forebrain_workspace$switchListForebrain,
    quantification = forebrain_workspace$salmonQuant
  ),
  Hindbrain = list(
    workspace = "foreBrain.RData",
    switch_object = forebrain_workspace$switchListHindbrain,
    quantification = NULL
  ),
  Midbrain = list(
    workspace = "isoform_final.RData",
    switch_object = midbrain_workspace$switchListMidbrain,
    quantification = midbrain_workspace$salmonQuant
  )
)

audit_rows <- lapply(names(regions), function(region) {
  item <- regions[[region]]
  switch_object <- item$switch_object
  retained_ids <- switch_object$isoformCountMatrix$isoform_id
  retained_profile <- refseq_profile(retained_ids)
  quantification_survives <- !is.null(item$quantification)

  if (quantification_survives) {
    quantification_ids <- item$quantification$abundance$isoform_id
    quantification_profile <- refseq_profile(quantification_ids)
    quantification_rows <- nrow(item$quantification$abundance)
    quantification_base_accessions <- quantification_profile$base_accessions
    quantification_suffix_instances <- quantification_profile$suffix_instances
    quantification_id_profile <- if (quantification_profile$all_match) {
      "All identifiers match NM_/NR_ RefSeq accession pattern with optional numeric suffix"
    } else {
      "At least one identifier falls outside the NM_/NR_ RefSeq accession pattern"
    }
    archived_import_options <- import_options(item$quantification)
    lineage <- if (count_lineage(item$quantification, switch_object)) {
      "Filtered salmonQuant counts exactly equal the retained isoformCountMatrix"
    } else {
      "No exact count-matrix lineage established"
    }
  } else {
    quantification_rows <- NA_integer_
    quantification_base_accessions <- NA_integer_
    quantification_suffix_instances <- NA_integer_
    quantification_id_profile <- "Separate pre-import salmonQuant object does not survive"
    archived_import_options <- NA_character_
    lineage <- "Only the post-import switchAnalyzeRlist count matrix survives"
  }

  exon_ids <- S4Vectors::mcols(switch_object$exons)$isoform_id
  data.frame(
    region = region,
    workspace = item$workspace,
    archived_input_stage = if (quantification_survives) {
      "salmonQuant plus post-import switchAnalyzeRlist"
    } else {
      "post-import switchAnalyzeRlist only"
    },
    salmonQuant_object_survives = quantification_survives,
    salmonQuant_transcript_entries = quantification_rows,
    salmonQuant_base_refseq_accessions = quantification_base_accessions,
    salmonQuant_numeric_suffix_instances = quantification_suffix_instances,
    salmonQuant_identifier_profile = quantification_id_profile,
    retained_count_matrix_rows = nrow(switch_object$isoformCountMatrix),
    retained_base_refseq_accessions = retained_profile$base_accessions,
    retained_numeric_suffix_instances = retained_profile$suffix_instances,
    retained_identifiers_all_nm_nr_pattern = retained_profile$all_match,
    retained_exon_records = length(switch_object$exons),
    retained_exon_transcripts = length(unique(exon_ids)),
    retained_seqlevel_count = length(GenomeInfoDb::seqlevels(switch_object$exons)),
    count_matrix_lineage = lineage,
    archived_import_options = archived_import_options,
    IsoformSwitchAnalyzeR_version = switch_object$runInfo$IsoformSwitchAnalyzeR$version,
    annotation_evidence = paste(
      "UCSC-style chromosome names and NM_/NR_ accessions;",
      "numeric suffixes occur for repeated accession instances;",
      "no MSTRG-style identifier survives"
    ),
    cross_workspace_reference_evidence = if (region %in% c("Forebrain", "Midbrain")) {
      cross_workspace_evidence
    } else {
      NA_character_
    },
    unresolved = paste(
      "Exact GTF/GFF release and checksum, Salmon executable and index,",
      "quant.sf paths, and StringTie2-to-Salmon hand-off are absent"
    ),
    check.names = FALSE
  )
})

audit <- do.call(rbind, audit_rows)
output_file <- file.path(paper_dir, "tables", "archived_input_structure_audit.csv")
write.csv(audit, output_file, row.names = FALSE, na = "")
message("Wrote ", normalizePath(output_file, winslash = "/", mustWork = TRUE))
