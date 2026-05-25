#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: 03.from_seurat_to_h5ad.r <input_seurat.rds> <output.h5ad>")
}

input_file <- args[1]
output_file <- args[2]

suppressPackageStartupMessages({
  library(Seurat)
  library(sceasy)
  library(reticulate)
})

# Optional: set RETICULATE_CONDA_ENV if a specific Python environment is needed.
reticulate_env <- Sys.getenv("RETICULATE_CONDA_ENV", unset = "")
if (nzchar(reticulate_env)) {
  reticulate::use_condaenv(reticulate_env, required = TRUE)
}

seurat_obj <- readRDS(input_file)
sceasy::convertFormat(
  seurat_obj,
  from = "seurat",
  to = "anndata",
  outFile = output_file
)
