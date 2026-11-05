#!/usr/bin/env Rscript
# Adapted from POTTR (https://github.com/AlBi-HHU/POTTR),
# data/TRACERx preprocessing/create_tree_files.r, MIT License, (c) 2025 Sara Schulte.
# Changes: reads the release files from one flat directory and drops the
# tables this step does not use.
#
# Writes every candidate tree reported by CONIPHER for each TRACERx patient
# to <output_path>/<patient>.xlsx, one sheet per tree. Requires the R package
# openxlsx.
#
# Usage: Rscript create_tree_files.r <input_path> <output_path>
#   <input_path> must contain 20221109_TRACERx421_phylogenetic_trees.rds.
library(openxlsx)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript create_tree_files.r <input_path> <output_path>")
}

input_path <- args[1]
outdir     <- args[2]

if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

data <- readRDS(file.path(input_path, "20221109_TRACERx421_phylogenetic_trees.rds"))

keys <- names(data)

for (i in seq_len(length(keys))) {
  print(i)
  key <- keys[i]
  # alt_trees holds all trees reported by CONIPHER; to keep only the default
  # tree, write data[[key]]$graph_pyclone$Corrected_tree instead.
  file_name <- file.path(outdir, paste0(key, '.xlsx'))
  write.xlsx(data[[key]]$graph_pyclone$alt_trees, file_name)
}
