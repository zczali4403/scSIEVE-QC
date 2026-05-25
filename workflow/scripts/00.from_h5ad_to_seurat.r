#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: 00.from_h5ad_to_seurat.r <input.h5ad> <output_seurat.rds>")
}

h5ad_file <- args[1]
seurat_file <- args[2]

suppressPackageStartupMessages({
  library(sceasy)
  library(reticulate)
})

# Optional: set RETICULATE_CONDA_ENV if a specific Python environment is needed.
reticulate_env <- Sys.getenv("RETICULATE_CONDA_ENV", unset = "")
if (nzchar(reticulate_env)) {
  reticulate::use_condaenv(reticulate_env, required = TRUE)
}

sceasy::convertFormat(
  h5ad_file,
  from = "anndata",
  to = "seurat",
  outFile = seurat_file
)
