# =============================================================================
# Phase 1 — Differential-expression analysis of GSE255075
# Script: 01_GSE255075_differential_expression.R
# =============================================================================
#
# Project
# -------
# Placental ACE2 and glucose transporter expression in gestational diabetes
# mellitus
#
# Author
# ------
# Md. Masudul Haque
#
# Purpose
# -------
# Reconstruct and analyze the GSE255075 placental RNA-seq dataset comparing
# normoglycemic control and gestational diabetes mellitus placentas.
#
# This phase performs count-matrix validation, DESeq2 differential-expression
# analysis, variance-stabilizing transformation, PCA quality control, and
# volcano-plot generation.
#
# Reproducibility boundary
# ------------------------
# This analysis uses the public GSE255075 transcriptomic dataset and does not
# contain or depend on the institutional patient-level placental cohort.
#
# Differential-expression inference is performed with DESeq2 and
# Benjamini-Hochberg multiple-testing correction. VST values are used only for
# visualization and downstream correlation analyses and are not used for
# differential-expression testing.
#
# Inputs
# ------
# data/raw/public/GSE255075/GSE255075_2.csv
#
# Outputs
# -------
# data/processed/GSE255075/dds_GSE255075.rds
# data/processed/GSE255075/vst_GSE255075.rds
# data/processed/GSE255075/sample_info_GSE255075.rds
# outputs/tables/GSE255075/DESeq2_GDM_vs_Normal_shrunk.csv
# figures/exploratory/GSE255075/GSE255075_PCA.png
# figures/exploratory/GSE255075/GSE255075_volcano_GOI.png
#
# Analysis-critical settings
# --------------------------
# Control samples: Normal1, Normal2, Normal3
# GDM samples: GDM1, GDM2, GDM3
# Genes retained when total count >= 10 across samples
# DESeq2 design: ~ condition
# FDR method: Benjamini-Hochberg
# Significance threshold: adjusted P < 0.05
# Random seed: 123
# =============================================================================

if (!requireNamespace("here", quietly = TRUE)) {
  stop("Install package 'here' before running this script.")
}
source(here::here("R", "00_project_setup.R"))

input_csv_255 <- file.path(path_data_raw_255, "GSE255075_2.csv")
assert_file_exists(input_csv_255)

gene_col_255 <- "GeneID"
control_samples <- c("Normal1", "Normal2", "Normal3")
gdm_samples <- c("GDM1", "GDM2", "GDM3")
analysis_samples <- c(control_samples, gdm_samples)

message("\n===== GSE255075: differential expression =====\n")

# -----------------------------------------------------------------------------
# 1. Load and validate count data
# -----------------------------------------------------------------------------

raw_data <- read.csv(
  input_csv_255,
  header = TRUE,
  sep = ",",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (any(colnames(raw_data) == "" | is.na(colnames(raw_data)))) {
  raw_data <- raw_data[
    ,
    !(colnames(raw_data) == "" | is.na(colnames(raw_data))),
    drop = FALSE
  ]
}

stopifnot(gene_col_255 %in% colnames(raw_data))
stopifnot(all(analysis_samples %in% colnames(raw_data)))

raw_data <- raw_data[
  ,
  c(analysis_samples, gene_col_255),
  drop = FALSE
]

raw_data <- raw_data[
  !is.na(raw_data[[gene_col_255]]) & raw_data[[gene_col_255]] != "",
  ,
  drop = FALSE
]

if (any(grepl("^ENSG", raw_data[[gene_col_255]]))) {
  raw_data[[gene_col_255]] <- sub(
    "\\..*$",
    "",
    raw_data[[gene_col_255]]
  )
}

raw_data <- raw_data |>
  dplyr::distinct(.data[[gene_col_255]], .keep_all = TRUE) |>
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(analysis_samples),
      as.numeric
    )
  )

count_vals <- as.matrix(raw_data[, analysis_samples, drop = FALSE])

if (anyNA(count_vals)) {
  stop("NA values were introduced in GSE255075 count columns.")
}

if (any(count_vals < 0)) {
  stop("Negative counts detected. DESeq2 requires non-negative raw counts.")
}

non_integer_fraction <- mean(
  abs(count_vals - round(count_vals)) > 1e-6
)

if (non_integer_fraction > 0.05) {
  warning(
    ">5% of values are non-integer. Confirm that the file contains raw counts."
  )
}

count_df_255 <- raw_data |>
  tibble::column_to_rownames(var = gene_col_255)

count_mat_255 <- as.matrix(count_df_255)
count_mat_255 <- count_mat_255[, analysis_samples, drop = FALSE]

if (any(abs(count_mat_255 - round(count_mat_255)) > 1e-6)) {
  warning("Non-integer values detected; rounding before DESeq2.")
}

count_mat_255 <- round(count_mat_255)
storage.mode(count_mat_255) <- "integer"

# -----------------------------------------------------------------------------
# 2. Sample metadata
# -----------------------------------------------------------------------------

sample_info_255 <- data.frame(
  row.names = analysis_samples,
  condition = factor(
    c(
      rep("Normal", length(control_samples)),
      rep("GDM", length(gdm_samples))
    ),
    levels = c("Normal", "GDM")
  )
)

# -----------------------------------------------------------------------------
# 3. DESeq2
# -----------------------------------------------------------------------------

dds_255 <- DESeq2::DESeqDataSetFromMatrix(
  countData = count_mat_255,
  colData = sample_info_255,
  design = ~ condition
)

dds_255 <- dds_255[
  rowSums(DESeq2::counts(dds_255)) >= 10,
]

dds_255 <- DESeq2::DESeq(dds_255)

res_255 <- DESeq2::results(
  dds_255,
  contrast = c("condition", "GDM", "Normal"),
  alpha = 0.05
)

coef_candidates <- grep(
  "^condition_.*GDM.*vs.*Normal$",
  DESeq2::resultsNames(dds_255),
  value = TRUE
)

if (length(coef_candidates) != 1) {
  stop(
    "Could not identify a unique GDM-vs-Normal coefficient for apeglm shrinkage.\n",
    "Available coefficients: ",
    paste(DESeq2::resultsNames(dds_255), collapse = ", ")
  )
}

res_shrunk_255 <- DESeq2::lfcShrink(
  dds_255,
  coef = coef_candidates,
  type = "apeglm"
)

res_df_255 <- as.data.frame(res_shrunk_255) |>
  tibble::rownames_to_column("gene") |>
  dplyr::arrange(padj)

write.csv(
  res_df_255,
  file.path(
    path_tables_255,
    "DESeq2_GDM_vs_Normal_shrunk.csv"
  ),
  row.names = FALSE
)

# Keep the unshrunk result as well because it is useful for diagnostics.
write.csv(
  as.data.frame(res_255) |>
    tibble::rownames_to_column("gene"),
  file.path(
    path_tables_255,
    "DESeq2_GDM_vs_Normal_unshrunk.csv"
  ),
  row.names = FALSE
)

deg_summary <- data.frame(
  criterion = c(
    "padj < 0.05",
    "padj < 0.05 and abs(log2FC) > 1"
  ),
  n = c(
    sum(!is.na(res_df_255$padj) & res_df_255$padj < 0.05),
    sum(
      !is.na(res_df_255$padj) &
        res_df_255$padj < 0.05 &
        abs(res_df_255$log2FoldChange) > 1
    )
  )
)

write.csv(
  deg_summary,
  file.path(path_tables_255, "DEG_summary.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 4. VST and persisted intermediate objects
# -----------------------------------------------------------------------------

vsd_255 <- DESeq2::vst(dds_255, blind = FALSE)

saveRDS(
  dds_255,
  file.path(path_processed_255, "dds_GSE255075.rds")
)

saveRDS(
  vsd_255,
  file.path(path_processed_255, "vst_GSE255075.rds")
)

saveRDS(
  sample_info_255,
  file.path(path_processed_255, "sample_info_GSE255075.rds")
)

# -----------------------------------------------------------------------------
# 5. PCA
# -----------------------------------------------------------------------------

pca_df_255 <- DESeq2::plotPCA(
  vsd_255,
  intgroup = "condition",
  returnData = TRUE
)

pct_var_255 <- round(
  100 * attr(pca_df_255, "percentVar")
)

p_pca_255 <- ggplot2::ggplot(
  pca_df_255,
  ggplot2::aes(PC1, PC2, color = condition)
) +
  ggplot2::geom_point(size = 3) +
  ggplot2::theme_bw() +
  ggplot2::labs(
    title = "PCA — GSE255075",
    x = paste0("PC1: ", pct_var_255[1], "%"),
    y = paste0("PC2: ", pct_var_255[2], "%"),
    color = "Condition"
  )

save_plot(
  p_pca_255,
  file.path(path_fig_255, "GSE255075_PCA.png"),
  width = 6.5,
  height = 5
)

# -----------------------------------------------------------------------------
# 6. Volcano plot
# -----------------------------------------------------------------------------

p_volcano_255 <- EnhancedVolcano::EnhancedVolcano(
  res_df_255,
  lab = res_df_255$gene,
  x = "log2FoldChange",
  y = "padj",
  selectLab = intersect(
    genes_of_interest,
    res_df_255$gene
  ),
  pCutoff = 0.05,
  FCcutoff = 1,
  xlim = c(-5, 5),
  ylim = c(-1, 30),
  pointSize = 2.0,
  labSize = 3.5,
  title = "GDM vs Normal — GSE255075"
)

ggplot2::ggsave(
  filename = file.path(
    path_fig_255,
    "GSE255075_volcano_GOI.png"
  ),
  plot = p_volcano_255,
  width = 8,
  height = 6.5,
  units = "in",
  dpi = dpi_out,
  bg = "white"
)

# -----------------------------------------------------------------------------
# 6B. Manuscript-style full volcano plot
# -----------------------------------------------------------------------------
#
# Reproduces the historical Figure 1B configuration.
# Uses the unshrunk DESeq2 result for the global transcriptomic volcano.
# Differential-expression inference remains BH-adjusted DESeq2 inference.
#
# Historical plot settings:
#   x = unshrunk DESeq2 log2 fold change
#   y = adjusted P value (padj)
#   FC cutoff = 1
#   adjusted-P cutoff = 1e-5
#   x range = -5 to 5
#   y range = 0 to 30
#

p_volcano_manuscript_255 <- EnhancedVolcano::EnhancedVolcano(
  res_255,
  lab = rownames(res_255),
  x = "log2FoldChange",
  y = "padj",

  pCutoff = 1e-5,
  FCcutoff = 1,

  xlim = c(-5, 5),
  ylim = c(0, 30),

  border = "full",
  borderWidth = 1.5,
  borderColour = "black",

  gridlines.major = FALSE,
  gridlines.minor = FALSE,

  title = "Control versus GDM",
  subtitle = NULL,

  legendLabels = c(
    "NS",
    expression(Log[2]~FC),
    "Adjusted P",
    expression("Adjusted P and"~Log[2]~FC)
  ),

  pointSize = 2,
  labSize = 3.5,
  max.overlaps = 20
)

ggplot2::ggsave(
  filename = file.path(
    path_fig_255,
    "GSE255075_volcano_manuscript.png"
  ),
  plot = p_volcano_manuscript_255,
  width = 6.7,
  height = 6.7,
  units = "in",
  dpi = dpi_out,
  bg = "white"
)
