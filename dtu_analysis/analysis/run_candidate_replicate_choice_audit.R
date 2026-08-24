#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 160)

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
    "Could not find both required .RData inputs. Pass their directory as the sole ",
    "command-line argument."
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
stages <- c("14.5", "15.5", "16.5")

extract_fraction_matrix <- function(x) {
  value <- x$isoformRepIF
  stopifnot(!anyDuplicated(value$isoform_id))
  answer <- as.matrix(value[, -1, drop = FALSE])
  storage.mode(answer) <- "numeric"
  rownames(answer) <- value$isoform_id
  answer
}
fractions <- lapply(objects, extract_fraction_matrix)

panel <- read.csv(
  file.path(table_dir, "preregistered_validation_panel.csv"),
  check.names = FALSE
)
required_accessions <- unique(c(
  panel$higher_midbrain_accession,
  panel$lower_midbrain_accession
))
for (region in regions) {
  absent <- setdiff(required_accessions, rownames(fractions[[region]]))
  if (length(absent)) {
    stop("Panel accessions absent from ", region, ": ", paste(absent, collapse = ", "))
  }
}

cell_keys <- as.vector(outer(regions, stages, paste, sep = "__"))
choice_grid <- expand.grid(rep(list(c("1", "2")), length(cell_keys)))
names(choice_grid) <- cell_keys
stopifnot(nrow(choice_grid) == 2^9)

cell_values <- function(accession, transform = identity) {
  answer <- matrix(
    NA_real_,
    nrow = length(regions),
    ncol = length(stages),
    dimnames = list(regions, stages)
  )
  for (region in regions) {
    for (stage in stages) {
      columns <- grep(
        paste0("^", gsub("\\.", "\\\\.", stage), "_[12]$"),
        colnames(fractions[[region]]),
        value = TRUE
      )
      stopifnot(length(columns) == 2L)
      answer[region, stage] <- mean(transform(fractions[[region]][accession, columns]))
    }
  }
  answer
}

choice_contrasts <- function(accession, transform = identity) {
  answer <- numeric(nrow(choice_grid))
  for (i in seq_len(nrow(choice_grid))) {
    curves <- setNames(numeric(length(regions)), regions)
    for (region in regions) {
      selected <- vapply(stages, function(stage) {
        replicate <- choice_grid[i, paste(region, stage, sep = "__")]
        column <- paste(stage, replicate, sep = "_")
        transform(fractions[[region]][accession, column])
      }, numeric(1))
      curves[[region]] <- selected[["15.5"]] -
        0.5 * (selected[["14.5"]] + selected[["16.5"]])
    }
    answer[[i]] <- curves[["Midbrain"]] -
      0.5 * (curves[["Forebrain"]] + curves[["Hindbrain"]])
  }
  answer
}

regional_curvature <- function(values) {
  curves <- values[, "15.5"] - 0.5 * (values[, "14.5"] + values[, "16.5"])
  unname(curves[["Midbrain"]] - 0.5 * (curves[["Forebrain"]] + curves[["Hindbrain"]]))
}

epsilon <- 0.005
logit <- function(value) qlogis(pmin(pmax(value, epsilon), 1 - epsilon))

results <- do.call(rbind, lapply(seq_len(nrow(panel)), function(i) {
  higher <- panel$higher_midbrain_accession[[i]]
  lower <- panel$lower_midbrain_accession[[i]]
  higher_logit <- choice_contrasts(higher, logit)
  lower_logit <- choice_contrasts(lower, logit)
  higher_raw <- choice_contrasts(higher)
  lower_raw <- choice_contrasts(lower)
  data.frame(
    gene = panel$gene[[i]],
    higher_accession = higher,
    lower_accession = lower,
    replicate_choice_combinations = nrow(choice_grid),
    higher_full_logit_regional_curvature = regional_curvature(cell_values(higher, logit)),
    higher_logit_choice_min = min(higher_logit),
    higher_logit_choice_max = max(higher_logit),
    higher_logit_fraction_positive = mean(higher_logit > 0),
    lower_full_logit_regional_curvature = regional_curvature(cell_values(lower, logit)),
    lower_logit_choice_min = min(lower_logit),
    lower_logit_choice_max = max(lower_logit),
    lower_logit_fraction_negative = mean(lower_logit < 0),
    joint_logit_expected_direction_fraction = mean(higher_logit > 0 & lower_logit < 0),
    higher_full_raw_regional_curvature = regional_curvature(cell_values(higher)),
    higher_raw_fraction_positive = mean(higher_raw > 0),
    lower_full_raw_regional_curvature = regional_curvature(cell_values(lower)),
    lower_raw_fraction_negative = mean(lower_raw < 0),
    joint_raw_expected_direction_fraction = mean(higher_raw > 0 & lower_raw < 0),
    stringsAsFactors = FALSE
  )
}))

stopifnot(
  nrow(results) == 6L,
  all(results$replicate_choice_combinations == 512L),
  all(is.finite(as.matrix(results[, 5:ncol(results)])))
)

write.csv(
  results,
  file.path(table_dir, "candidate_replicate_choice_audit.csv"),
  row.names = FALSE
)

cat("Candidate replicate-choice audit complete.\n")
print(results)
