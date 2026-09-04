# =============================================================================
# Phase 0 — Shared project configuration and helper functions
# Script: 00_project_setup.R
# =============================================================================
#
# Project
# -------
# Placental ACE2 and glucose transporter expression in gestational diabetes
# mellitus
#
# Purpose
# -------
# Define project-relative paths, required packages, shared gene panels,
# plotting defaults, reproducibility settings, and helper functions used by
# all downstream analysis scripts.
#
# Reproducibility boundary
# ------------------------
# This script performs no scientific hypothesis testing. It establishes a
# common execution environment so that all downstream scripts can be run from
# the project root on macOS, Windows, or Linux without machine-specific paths.
#
# Inputs
# ------
# None.
#
# Outputs
# -------
# Creates the local processed-data, output, log, and figure directories when
# they do not already exist.
#
# Analysis-critical settings
# --------------------------
# Random seed: 123
# Default raster figure resolution: 600 dpi
# =============================================================================

if (!requireNamespace("here", quietly = TRUE)) {
  stop(
    "Package 'here' is required. Install it with install.packages('here') ",
    "and rerun the script."
  )
}

required_packages <- c(
  "readxl",
  "DESeq2",
  "apeglm",
  "EnhancedVolcano",
  "gage",
  "pathview",
  "clusterProfiler",
  "org.Hs.eg.db",
  "enrichplot",
  "pheatmap",
  "ComplexHeatmap",
  "circlize",
  "ggplot2",
  "dplyr",
  "tidyr",
  "tibble",
  "stringr",
  "AnnotationDbi",
  "openxlsx"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required R packages: ",
    paste(missing_packages, collapse = ", "),
    "\nInstall them before running the analysis."
  )
}

suppressPackageStartupMessages({
  library(here)
  library(readxl)
  library(DESeq2)
  library(apeglm)
  library(EnhancedVolcano)
  library(gage)
  library(pathview)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(pheatmap)
  library(ComplexHeatmap)
  library(circlize)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(AnnotationDbi)
  library(openxlsx)
})

set.seed(123)

# -----------------------------------------------------------------------------
# Project paths
# -----------------------------------------------------------------------------

path_data_raw_255 <- here("data", "raw", "public", "GSE255075")
path_data_raw_249 <- here("data", "raw", "public", "GSE249311")
path_data_private <- here("data", "private", "clinical")

path_processed_255 <- here("data", "processed", "GSE255075")
path_processed_249 <- here("data", "processed", "GSE249311")

path_tables_255 <- here("outputs", "tables", "GSE255075")
path_tables_249 <- here("outputs", "tables", "GSE249311")
path_stats_clinical <- here("outputs", "statistics", "clinical")
path_logs <- here("outputs", "logs")

path_fig_main <- here("figures", "main")
path_fig_supp <- here("figures", "supplementary")
path_fig_255 <- here("figures", "exploratory", "GSE255075")
path_fig_249 <- here("figures", "exploratory", "GSE249311")

directories_to_create <- c(
  path_processed_255,
  path_processed_249,
  path_tables_255,
  path_tables_249,
  path_stats_clinical,
  path_logs,
  path_fig_main,
  path_fig_supp,
  path_fig_255,
  path_fig_249
)

invisible(lapply(
  directories_to_create,
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

# -----------------------------------------------------------------------------
# Shared analysis configuration
# -----------------------------------------------------------------------------

dpi_out <- 600
box_line <- 0.8
point_size <- 3

glut_genes <- c("SLC2A1", "SLC2A3", "SLC2A4", "SLC2A9")

ras_core <- c(
  "ACE", "ACE2", "REN", "AGT", "AGTR1", "AGTR2", "MAS1"
)

ras_reg <- c(
  "MME", "THOP1", "ANPEP", "LNPEP", "ADAM17", "TIMP3"
)

ras_genes <- unique(c(ras_core, ras_reg))
gene_panel <- unique(c(glut_genes, ras_core, ras_reg))
genes_of_interest <- unique(c(glut_genes, ras_genes))

# Correlation panel retained from the historical analysis.
reg_genes_cor <- c(
  "REN", "AGT", "ACE", "ACE2", "AGTR1", "MME", "THOP1",
  "ANPEP", "LNPEP", "ADAM17", "TIMP3", "MAS1"
)

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

save_plot <- function(
    p,
    filepath,
    width = 8,
    height = 5,
    dpi = dpi_out,
    bg = "white"
) {
  dir.create(dirname(filepath), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    filename = filepath,
    plot = p,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    bg = bg
  )
}

safe_name <- function(x, maxlen = 80) {
  x <- stringr::str_replace_all(x, "[^A-Za-z0-9_\\-]+", "_")
  x <- stringr::str_replace_all(x, "_+", "_")
  x <- stringr::str_replace_all(x, "^_|_$", "")
  substr(x, 1, maxlen)
}

write_session_info <- function(script_name) {
  outfile <- file.path(
    path_logs,
    paste0(tools::file_path_sans_ext(basename(script_name)), "_sessionInfo.txt")
  )
  capture.output(sessionInfo(), file = outfile)
  invisible(outfile)
}

assert_file_exists <- function(path) {
  if (!file.exists(path)) {
    stop(
      "Required input file not found:\n  ", path,
      "\nSee README.md and data/README.md for the expected project structure."
    )
  }
  invisible(path)
}

assert_object_exists <- function(path, description = "processed object") {
  if (!file.exists(path)) {
    stop(
      "Required ", description, " not found:\n  ", path,
      "\nRun the preceding analysis phase first."
    )
  }
  invisible(path)
}

message("Project root: ", here())
