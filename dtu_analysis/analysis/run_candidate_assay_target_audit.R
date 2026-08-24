#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 160)

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(Biostrings)
  library(BSgenome.Mmusculus.UCSC.mm10)
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

required_inputs <- "isoform_final.RData"
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) {
  stop("Expected at most one argument: the directory containing isoform_final.RData.")
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
    "Could not find isoform_final.RData. Pass its directory as the sole ",
    "command-line argument."
  )
}
input_dir <- normalizePath(
  input_candidates[which(complete_input)[[1]]],
  winslash = "/",
  mustWork = TRUE
)
message("Loading archived input from: ", input_dir)
workspace <- new.env(parent = globalenv())
load(file.path(input_dir, "isoform_final.RData"), envir = workspace)
if (!exists("switchListMidbrain", envir = workspace, inherits = FALSE)) {
  stop("switchListMidbrain is absent from isoform_final.RData.")
}

top_candidates <- read.csv(
  file.path(table_dir, "transient_regional_top_candidates.csv"),
  check.names = FALSE
)
top_isoforms <- read.csv(
  file.path(table_dir, "transient_regional_top_candidate_isoforms.csv"),
  check.names = FALSE
)
panel <- data.frame(
  gene = top_candidates$gene_id,
  calibrated_scan_rank = seq_len(nrow(top_candidates)),
  higher_midbrain_accession = vapply(top_candidates$gene_id, function(gene) {
    value <- top_isoforms$isoform_id[
      top_isoforms$gene_id == gene & top_isoforms$direction == "higher"
    ]
    if (length(value) != 1L) stop("Expected one higher accession for ", gene)
    value
  }, character(1)),
  lower_midbrain_accession = vapply(top_candidates$gene_id, function(gene) {
    value <- top_isoforms$isoform_id[
      top_isoforms$gene_id == gene & top_isoforms$direction == "lower"
    ]
    if (length(value) != 1L) stop("Expected one lower accession for ", gene)
    value
  }, character(1)),
  stringsAsFactors = FALSE
)
stopifnot(
  nrow(panel) == 6L,
  !anyDuplicated(panel$gene),
  !anyNA(panel$higher_midbrain_accession),
  !anyNA(panel$lower_midbrain_accession)
)

exons <- as.data.frame(workspace$switchListMidbrain$exons)
required_columns <- c(
  "seqnames", "start", "end", "strand", "isoform_id", "gene_id", "gene_name"
)
if (!all(required_columns %in% colnames(exons))) {
  stop("The archived exon model lacks required columns: ", paste(
    setdiff(required_columns, colnames(exons)), collapse = ", "
  ))
}
exons <- unique(exons[
  exons$gene_id %in% panel$gene,
  required_columns,
  drop = FALSE
])

panel_long <- rbind(
  data.frame(
    gene = panel$gene,
    calibrated_scan_rank = panel$calibrated_scan_rank,
    expected_direction = "higher",
    accession = panel$higher_midbrain_accession,
    paired_accession = panel$lower_midbrain_accession
  ),
  data.frame(
    gene = panel$gene,
    calibrated_scan_rank = panel$calibrated_scan_rank,
    expected_direction = "lower",
    accession = panel$lower_midbrain_accession,
    paired_accession = panel$higher_midbrain_accession
  )
)
panel_long$expected_direction <- factor(
  panel_long$expected_direction,
  levels = c("higher", "lower")
)
panel_long <- panel_long[order(
  panel_long$calibrated_scan_rank,
  panel_long$expected_direction
), ]
panel_long$expected_direction <- as.character(panel_long$expected_direction)

missing_accessions <- setdiff(panel_long$accession, exons$isoform_id)
if (length(missing_accessions)) {
  stop("Panel accessions absent from the archived exon model: ", paste(
    missing_accessions, collapse = ", "
  ))
}

ordered_exons <- function(accession) {
  value <- exons[exons$isoform_id == accession, , drop = FALSE]
  value <- value[order(value$start, value$end), , drop = FALSE]
  if (length(unique(value$seqnames)) != 1L ||
      length(unique(value$strand)) != 1L ||
      length(unique(value$gene_id)) != 1L) {
    stop("Inconsistent locus metadata for ", accession)
  }
  if (any(value$start[-1] <= head(value$end, -1))) {
    stop("Overlapping or unsorted exons prevent junction construction for ", accession)
  }
  n_exons <- nrow(value)
  value$genomic_exon_rank <- seq_len(n_exons)
  value$transcript_exon_rank_5to3 <- if (value$strand[[1]] == "+") {
    seq_len(n_exons)
  } else {
    rev(seq_len(n_exons))
  }
  value
}

reference_genome <- BSgenome.Mmusculus.UCSC.mm10
gene_isoforms <- unique(exons$isoform_id)
gene_by_isoform <- setNames(
  exons$gene_id[match(gene_isoforms, exons$isoform_id)],
  gene_isoforms
)

transcript_exon_sequences <- setNames(lapply(gene_isoforms, function(accession) {
  value <- ordered_exons(accession)
  if (value$strand[[1]] == "-") {
    value <- value[rev(seq_len(nrow(value))), , drop = FALSE]
  }
  ranges <- GRanges(
    seqnames = value$seqnames,
    ranges = IRanges(value$start, value$end),
    strand = value$strand
  )
  as.character(getSeq(reference_genome, ranges))
}), gene_isoforms)
transcript_sequences <- vapply(
  transcript_exon_sequences,
  paste0,
  collapse = "",
  FUN.VALUE = character(1)
)

archived_sequences <- workspace$combinedMitbrain$ntSequence
sequence_audit <- do.call(rbind, lapply(seq_len(nrow(panel_long)), function(i) {
  item <- panel_long[i, ]
  accession <- item$accession
  reconstructed <- transcript_sequences[[accession]]
  survives <- accession %in% names(archived_sequences)
  archived <- if (survives) as.character(archived_sequences[[accession]]) else NA_character_
  data.frame(
    gene = item$gene,
    calibrated_scan_rank = item$calibrated_scan_rank,
    expected_direction = item$expected_direction,
    accession = accession,
    paired_accession = item$paired_accession,
    reconstructed_transcript_length_nt = nchar(reconstructed),
    reconstructed_ambiguous_bases = lengths(regmatches(
      reconstructed,
      gregexpr("N", reconstructed, fixed = TRUE)
    )),
    archived_sequence_survives = survives,
    archived_sequence_length_nt = if (survives) nchar(archived) else NA_integer_,
    exact_match_to_archived_sequence = if (survives) identical(reconstructed, archived) else NA,
    sequence_evidence = if (survives) {
      "Exact exon-plus-mm10 reconstruction match to archived ntSequence"
    } else {
      "Reconstructed from archived exons and named mm10 reference; no candidate ntSequence survives"
    },
    reference_package = "BSgenome.Mmusculus.UCSC.mm10",
    reference_package_version = as.character(packageVersion(
      "BSgenome.Mmusculus.UCSC.mm10"
    )),
    reference_genome = metadata(reference_genome)$genome,
    annotation_boundary = paste(
      "Sequence reconstruction is auditable, but the original RefSeq annotation",
      "release and accession version suffix remain unresolved"
    ),
    stringsAsFactors = FALSE
  )
}))
stopifnot(
  nrow(sequence_audit) == 12L,
  sum(sequence_audit$archived_sequence_survives) == 10L,
  all(sequence_audit$exact_match_to_archived_sequence[
    sequence_audit$archived_sequence_survives
  ]),
  all(sequence_audit$reconstructed_ambiguous_bases == 0L)
)

make_junctions <- function(accession) {
  value <- ordered_exons(accession)
  n_exons <- nrow(value)
  if (n_exons < 2L) {
    return(data.frame())
  }
  index <- seq_len(n_exons - 1L)
  strand <- value$strand[[1]]
  left_end <- value$end[index]
  right_start <- value$start[index + 1L]
  transcript_rank <- if (strand == "+") index else n_exons - index
  data.frame(
    gene = value$gene_id[[1]],
    accession = accession,
    chromosome = as.character(value$seqnames[[1]]),
    strand = strand,
    archived_exon_count = n_exons,
    transcript_junction_rank_5to3 = transcript_rank,
    genomic_left_exon_start_1based = value$start[index],
    genomic_left_exon_end_1based = left_end,
    genomic_right_exon_start_1based = right_start,
    genomic_right_exon_end_1based = value$end[index + 1L],
    donor_exonic_base_1based = if (strand == "+") left_end else right_start,
    acceptor_exonic_base_1based = if (strand == "+") right_start else left_end,
    archived_intron_length_nt = right_start - left_end - 1L,
    junction_key = paste(
      value$seqnames[[1]], left_end, right_start, strand, sep = ":"
    ),
    stringsAsFactors = FALSE
  )
}

all_gene_junctions <- do.call(rbind, lapply(gene_isoforms, make_junctions))
panel_junctions <- do.call(rbind, lapply(panel_long$accession, make_junctions))
panel_junctions <- merge(
  panel_junctions,
  panel_long[, c(
    "gene", "calibrated_scan_rank", "expected_direction", "accession",
    "paired_accession"
  )],
  by = c("gene", "accession"),
  all.x = TRUE,
  sort = FALSE
)

panel_junctions$pair_discriminating <- mapply(
  function(key, paired_accession) {
    !key %in% all_gene_junctions$junction_key[
      all_gene_junctions$accession == paired_accession
    ]
  },
  panel_junctions$junction_key,
  panel_junctions$paired_accession
)
panel_junctions$archived_gene_isoforms_with_junction <- mapply(
  function(gene, key) {
    length(unique(all_gene_junctions$accession[
      all_gene_junctions$gene == gene & all_gene_junctions$junction_key == key
    ]))
  },
  panel_junctions$gene,
  panel_junctions$junction_key
)
panel_junctions$archived_gene_unique_junction <-
  panel_junctions$archived_gene_isoforms_with_junction == 1L
panel_junctions$preferred_coordinate_candidate <-
  panel_junctions$pair_discriminating &
  panel_junctions$archived_gene_unique_junction

count_fixed_occurrences <- function(pattern, subject) {
  locations <- gregexpr(pattern, subject, fixed = TRUE)[[1]]
  if (locations[[1]] == -1L) 0L else length(locations)
}
junction_kmers <- mapply(
  function(accession, junction_rank) {
    exon_sequences <- transcript_exon_sequences[[accession]]
    rank <- as.integer(junction_rank)
    if (nchar(exon_sequences[[rank]]) < 20L ||
        nchar(exon_sequences[[rank + 1L]]) < 20L) {
      stop("A 20+20 junction k-mer cannot be constructed for ", accession)
    }
    paste0(
      substr(
        exon_sequences[[rank]],
        nchar(exon_sequences[[rank]]) - 19L,
        nchar(exon_sequences[[rank]])
      ),
      substr(exon_sequences[[rank + 1L]], 1L, 20L)
    )
  },
  panel_junctions$accession,
  panel_junctions$transcript_junction_rank_5to3,
  USE.NAMES = FALSE
)
panel_junctions$junction_kmer_5to3_20plus20 <- junction_kmers
panel_junctions$junction_kmer_length_nt <- nchar(junction_kmers)
panel_junctions$junction_kmer_gc_fraction <- vapply(
  junction_kmers,
  function(value) {
    bases <- strsplit(value, "", fixed = TRUE)[[1]]
    mean(bases %in% c("G", "C"))
  },
  numeric(1)
)
panel_junctions$junction_kmer_occurrences_in_paired_accession <- mapply(
  function(kmer, paired_accession) {
    count_fixed_occurrences(kmer, transcript_sequences[[paired_accession]])
  },
  junction_kmers,
  panel_junctions$paired_accession
)
panel_junctions$junction_kmer_occurrences_in_archived_same_gene_models <- mapply(
  function(kmer, gene) {
    sum(vapply(
      transcript_sequences[gene_by_isoform == gene],
      function(subject) count_fixed_occurrences(kmer, subject),
      FUN.VALUE = integer(1)
    ))
  },
  junction_kmers,
  panel_junctions$gene
)
panel_junctions$sequence_candidate_status <- ifelse(
  panel_junctions$preferred_coordinate_candidate &
    panel_junctions$junction_kmer_occurrences_in_paired_accession == 0L &
    panel_junctions$junction_kmer_occurrences_in_archived_same_gene_models == 1L,
  "Coordinate- and 40-nt sequence-discriminating within archived same-gene models",
  ifelse(
    panel_junctions$preferred_coordinate_candidate,
    "Coordinate-discriminating but 40-nt sequence candidate is not unique within archived same-gene models",
    "Not a preferred discriminating junction"
  )
)
panel_junctions$assay_interpretation <- ifelse(
  panel_junctions$preferred_coordinate_candidate,
  "Archived-gene-unique splice junction; candidate for accession-discriminating assay design",
  ifelse(
    panel_junctions$pair_discriminating,
    "Discriminates the frozen pair but is shared by another archived isoform of the gene",
    "Shared by both frozen accessions; not a discriminating target"
  )
)
panel_junctions$coordinate_status <- paste(
  "Archived UCSC-style coordinates with missing genome/annotation version;",
  "revalidate against the frozen prospective reference before oligo design"
)
panel_junctions <- panel_junctions[order(
  panel_junctions$calibrated_scan_rank,
  match(panel_junctions$expected_direction, c("higher", "lower")),
  panel_junctions$transcript_junction_rank_5to3
), ]

make_unique_segments <- function(accession, gene) {
  target_exons <- ordered_exons(accession)
  target_ranges <- GRanges(
    seqnames = target_exons$seqnames,
    ranges = IRanges(target_exons$start, target_exons$end),
    strand = target_exons$strand
  )
  other <- exons[exons$gene_id == gene & exons$isoform_id != accession, ]
  other_ranges <- GRanges(
    seqnames = other$seqnames,
    ranges = IRanges(other$start, other$end),
    strand = other$strand
  )
  unique_ranges <- GenomicRanges::setdiff(
    GenomicRanges::reduce(target_ranges),
    GenomicRanges::reduce(other_ranges)
  )
  if (!length(unique_ranges)) {
    return(data.frame())
  }
  exon_ranges <- GRanges(
    seqnames = target_exons$seqnames,
    ranges = IRanges(target_exons$start, target_exons$end),
    strand = target_exons$strand
  )
  overlaps <- findOverlaps(unique_ranges, exon_ranges)
  exon_rank <- vapply(seq_along(unique_ranges), function(i) {
    hits <- subjectHits(overlaps)[queryHits(overlaps) == i]
    paste(sort(unique(target_exons$transcript_exon_rank_5to3[hits])), collapse = ";")
  }, character(1))
  rank_numeric <- suppressWarnings(as.integer(exon_rank))
  terminal_status <- ifelse(
    rank_numeric == 1L,
    "5prime_terminal_exon_segment",
    ifelse(
      rank_numeric == nrow(target_exons),
      "3prime_terminal_exon_segment",
      "internal_exon_segment"
    )
  )
  data.frame(
    gene = gene,
    accession = accession,
    chromosome = as.character(seqnames(unique_ranges)),
    strand = as.character(strand(unique_ranges)),
    start_1based = start(unique_ranges),
    end_1based = end(unique_ranges),
    width_nt = width(unique_ranges),
    transcript_exon_rank_5to3 = exon_rank,
    segment_location = terminal_status,
    segment_key = paste(
      as.character(seqnames(unique_ranges)), start(unique_ranges),
      end(unique_ranges), as.character(strand(unique_ranges)), sep = ":"
    ),
    coordinate_status = paste(
      "Archived-gene-unique exonic segment with missing genome/annotation version;",
      "revalidate sequence and specificity before oligo design"
    ),
    stringsAsFactors = FALSE
  )
}

unique_segments <- do.call(rbind, lapply(seq_len(nrow(panel_long)), function(i) {
  make_unique_segments(panel_long$accession[[i]], panel_long$gene[[i]])
}))
if (nrow(unique_segments)) {
  unique_segments <- merge(
    unique_segments,
    panel_long[, c(
      "gene", "calibrated_scan_rank", "expected_direction", "accession",
      "paired_accession"
    )],
    by = c("gene", "accession"),
    all.x = TRUE,
    sort = FALSE
  )
  unique_segments <- unique_segments[order(
    unique_segments$calibrated_scan_rank,
    match(unique_segments$expected_direction, c("higher", "lower")),
    unique_segments$start_1based
  ), ]
}

summary_rows <- do.call(rbind, lapply(seq_len(nrow(panel_long)), function(i) {
  item <- panel_long[i, ]
  junctions <- panel_junctions[
    panel_junctions$accession == item$accession,
    ,
    drop = FALSE
  ]
  segments <- unique_segments[
    unique_segments$accession == item$accession,
    ,
    drop = FALSE
  ]
  preferred_junctions <- sum(junctions$preferred_coordinate_candidate)
  sequence_supported_junctions <- sum(
    junctions$sequence_candidate_status ==
      "Coordinate- and 40-nt sequence-discriminating within archived same-gene models"
  )
  usable_segments <- sum(segments$width_nt >= 20L)
  sequence_row <- sequence_audit[sequence_audit$accession == item$accession, ]
  target_class <- if (preferred_junctions > 0L) {
    "archived-gene-unique splice junction"
  } else if (usable_segments > 0L) {
    "archived-gene-unique exonic segment fallback"
  } else {
    "no archived coordinate-discriminating target"
  }
  data.frame(
    gene = item$gene,
    calibrated_scan_rank = item$calibrated_scan_rank,
    expected_direction = item$expected_direction,
    accession = item$accession,
    paired_accession = item$paired_accession,
    archived_exons = if (nrow(junctions)) junctions$archived_exon_count[[1]] else 1L,
    archived_junctions = nrow(junctions),
    pair_discriminating_junctions = sum(junctions$pair_discriminating),
    archived_gene_unique_junctions = sum(junctions$archived_gene_unique_junction),
    preferred_junction_coordinate_candidates = preferred_junctions,
    preferred_junction_sequence_candidates = sequence_supported_junctions,
    reconstructed_transcript_length_nt = sequence_row$reconstructed_transcript_length_nt,
    archived_sequence_survives = sequence_row$archived_sequence_survives,
    exact_match_to_archived_sequence = sequence_row$exact_match_to_archived_sequence,
    archived_gene_unique_exonic_segments = nrow(segments),
    gene_unique_segments_at_least_20nt = usable_segments,
    longest_gene_unique_segment_nt = if (nrow(segments)) max(segments$width_nt) else 0L,
    recommended_archived_target_class = target_class,
    assay_design_status = paste(
      "Coordinate audit only; no primer or probe is frozen until annotation,",
      "reference sequence and in-silico specificity are prospectively approved"
    ),
    stringsAsFactors = FALSE
  )
}))

armc8_lower <- summary_rows[
  summary_rows$gene == "Armc8" & summary_rows$expected_direction == "lower",
]
stopifnot(
  nrow(summary_rows) == 12L,
  sum(summary_rows$preferred_junction_coordinate_candidates > 0L) == 11L,
  sum(summary_rows$preferred_junction_sequence_candidates > 0L) == 11L,
  all(summary_rows$preferred_junction_sequence_candidates ==
        summary_rows$preferred_junction_coordinate_candidates),
  nrow(armc8_lower) == 1L,
  armc8_lower$preferred_junction_coordinate_candidates == 0L,
  armc8_lower$gene_unique_segments_at_least_20nt > 0L,
  all(summary_rows$recommended_archived_target_class !=
        "no archived coordinate-discriminating target"),
  all(panel_junctions$junction_kmer_length_nt == 40L),
  all(panel_junctions$junction_kmer_occurrences_in_paired_accession[
    panel_junctions$preferred_coordinate_candidate
  ] == 0L),
  all(panel_junctions$junction_kmer_occurrences_in_archived_same_gene_models[
    panel_junctions$preferred_coordinate_candidate
  ] == 1L)
)

write.csv(
  panel_junctions,
  file.path(table_dir, "candidate_junction_target_audit.csv"),
  row.names = FALSE
)
write.csv(
  unique_segments,
  file.path(table_dir, "candidate_unique_exonic_segments.csv"),
  row.names = FALSE
)
write.csv(
  summary_rows,
  file.path(table_dir, "candidate_assay_target_summary.csv"),
  row.names = FALSE
)
write.csv(
  sequence_audit,
  file.path(table_dir, "candidate_sequence_reconstruction_audit.csv"),
  row.names = FALSE
)

cat("Candidate assay-target audit complete.\n")
print(summary_rows)
