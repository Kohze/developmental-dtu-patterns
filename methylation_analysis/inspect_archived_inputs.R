#!/usr/bin/env Rscript

# Read-only inventory used by the professor-style reconstruction audit.
suppressPackageStartupMessages(library(GenomicRanges))

file_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script <- sub("^--file=", "", file_argument[1])
root <- normalizePath(file.path(dirname(script), "..", ".."),
                      winslash = "/", mustWork = FALSE)
if (!file.exists(file.path(root, "methylation_subset_grl.RData"))) {
  root <- normalizePath(file.path(getwd(), "..", ".."), winslash = "/", mustWork = TRUE)
}

for (filename in c("forebrain_isoforms_data_with_dmr.rds",
                   "midbrain_isoforms_data_with_dmr.rds")) {
  object <- readRDS(file.path(root, filename))
  cat("\n===", filename, "===\n")
  cat("class:", paste(class(object), collapse = "/"), "\n")
  cat("items:", length(object), "\n")
  cat("regions:\n")
  print(table(vapply(object, function(item) item$region, character(1))))
  item <- object[[1]]
  cat("first item:", item$isoform_id, item$region, "\n")
  cat("range:", as.character(item$range), "\n")
  cat("extended range:", as.character(item$extended_range), "\n")
  cat("methylation entries:", length(item$methylation), "\n")
  cat("methylation names:", paste(names(item$methylation), collapse = ","), "\n")
}

environment <- new.env(parent = emptyenv())
load(file.path(root, "methylation_subset_grl.RData"), envir = environment)
samples <- environment$methylation_subset_grl
cat("\n=== methylation_subset_grl.RData ===\n")
cat("samples:", length(samples), "\n")
cat("sample IDs:", paste(names(samples), collapse = ","), "\n")
cat("CpG calls per sample:\n")
print(summary(lengths(samples)))

for (filename in c("foreBrain.RData", "isoform_final.RData")) {
  environment <- new.env()
  load(file.path(root, filename), envir = environment)
  object_name <- if (filename == "foreBrain.RData") "combinedForebrain" else "combinedMitbrain"
  object <- environment[[object_name]]
  cat("\n===", filename, object_name, "===\n")
  cat("components:", paste(names(object), collapse = ","), "\n")
  for (component in c("isoformRepIF", "isoformRepExpression", "exons")) {
    value <- object[[component]]
    cat(component, "class:", paste(class(value), collapse = "/"),
        "dim:", paste(dim(value), collapse = "x"), "\n")
    if (!is.null(colnames(value))) {
      cat(component, "columns:", paste(colnames(value), collapse = ","), "\n")
    }
    print(utils::head(value, 2))
  }
  rm(environment, object)
  invisible(gc())
}
