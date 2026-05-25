#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: 01.from_seurat_to_sce.r <input_seurat.rds> <output_sce.rds>")
}

seurat_file <- args[1]
sce_file <- args[2]

suppressPackageStartupMessages({
  library(sceasy)
  library(reticulate)
})

# Optional: set RETICULATE_CONDA_ENV if a specific Python environment is needed.
reticulate_env <- Sys.getenv("RETICULATE_CONDA_ENV", unset = "")
if (nzchar(reticulate_env)) {
  reticulate::use_condaenv(reticulate_env, required = TRUE)
}

seurat_object <- readRDS(seurat_file)
sceasy::convertFormat(
  seurat_object,
  from = "seurat",
  to = "sce",
  outFile = sce_file
)
