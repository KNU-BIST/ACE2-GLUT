# =============================================================================
# Phase 4 — Differential-expression analysis of GSE249311
# Script: 04_GSE249311_differential_expression.R
# =============================================================================
#
# Project
# -------
# Placental ACE2 and glucose transporter expression in gestational diabetes
# mellitus
#
# Purpose
# -------
# Process and analyze placental RNA-seq counts from GSE249311 to provide
# independent transcriptomic context for glucose-transporter and
# renin-angiotensin-system expression across maternal glycemic categories.
#
# Reproducibility boundary
# ------------------------
# Clinical categories are retained according to the source dataset:
#   - Control
#   - GDMA1: lifestyle-managed gestational diabetes
#   - GDMA2: medication-requiring gestational diabetes
#   - T2DM: pregestational type 2 diabetes
#
# GDMA2 is not interpreted as synonymous with insulin treatment.
# Aspirin-exposure information is unavailable in this public dataset.
#
# Inputs
# ------
# data/raw/public/GSE249311/GSE249311_count_matrix.xlsx
#
# Outputs
# -------
# data/processed/GSE249311/dds_GSE249311.rds
# data/processed/GSE249311/vst_GSE249311.rds
# data/processed/GSE249311/sample_info_GSE249311.rds
# outputs/tables/GSE249311/DESeq2_*.csv
#
# Analysis-critical settings
# --------------------------
# Sample sizes:
#   Control = 11
#   GDMA1   = 5
#   GDMA2   = 9
#   T2DM    = 5
#
# The condition vector assumes the columns in the validated count matrix occur
# in exactly that order.
#
# DESeq2 design: ~ condition
# Multiple-testing correction: Benjamini-Hochberg
# =============================================================================

if (!requireNamespace("here", quietly = TRUE)) {
  stop("Install package 'here' before running this script.")
}
source(here::here("R", "00_project_setup.R"))

input_xlsx_249 <- file.path(
  path_data_raw_249,
  "GSE249311_count_matrix.xlsx"
)

assert_file_exists(input_xlsx_249)

condition_249 <- c(
  rep("Control", 11),
  rep("GDMA1", 5),
  rep("GDMA2", 9),
  rep("T2DM", 5)
)

message("\n===== GSE249311: differential expression =====\n")

# -----------------------------------------------------------------------------
# 1. Load and validate counts
# -----------------------------------------------------------------------------

counts_df_249 <- readxl::read_excel(
  input_xlsx_249
)

if (ncol(counts_df_249) < 2) {
  stop(
    "GSE249311 count matrix must contain one gene column plus sample columns."
  )
}

gene_ids_249 <- trimws(
  as.character(counts_df_249[[1]])
)

valid_gene <- !is.na(gene_ids_249) &
  gene_ids_249 != ""

counts_df_249 <- counts_df_249[
  valid_gene,
  ,
  drop = FALSE
]

gene_ids_249 <- gene_ids_249[
  valid_gene
]

count_data_249 <- counts_df_249[
  ,
  -1,
  drop = FALSE
]

count_data_249[] <- lapply(
  count_data_249,
  function(x) {
    suppressWarnings(
      as.numeric(as.character(x))
    )
  }
)

if (anyNA(count_data_249)) {
  bad_columns <- names(count_data_249)[
    vapply(
      count_data_249,
      anyNA,
      logical(1)
    )
  ]

  stop(
    "Missing/non-numeric count values found in: ",
    paste(bad_columns, collapse = ", ")
  )
}

counts_raw_249 <- as.matrix(
  count_data_249
)

rownames(counts_raw_249) <- gene_ids_249

counts_249 <- rowsum(
  counts_raw_249,
  group = rownames(counts_raw_249),
  reorder = FALSE
)

if (any(counts_249 < 0)) {
  stop("Negative count values found.")
}

if (any(abs(counts_249 - round(counts_249)) > 1e-8)) {
  warning(
    "Non-integer count values detected; rounding before DESeq2."
  )
}

counts_249 <- round(counts_249)
storage.mode(counts_249) <- "integer"

sample_names_249 <- colnames(counts_249)

if (length(sample_names_249) != length(condition_249)) {
  stop(
    "Expected 30 GSE249311 sample columns but found ",
    length(sample_names_249),
    "."
  )
}

# -----------------------------------------------------------------------------
# 2. Sample metadata
# -----------------------------------------------------------------------------

sample_info_249 <- data.frame(
  row.names = sample_names_249,
  condition = factor(
    condition_249,
    levels = c(
      "Control",
      "GDMA1",
      "GDMA2",
      "T2DM"
    )
  )
)

print(
  table(sample_info_249$condition)
)

expected_counts_249 <- c(
  Control = 11,
  GDMA1 = 5,
  GDMA2 = 9,
  T2DM = 5
)

if (!all(
  table(sample_info_249$condition) ==
    expected_counts_249
)) {
  stop(
    "GSE249311 group counts do not match the expected 11/5/9/5 structure."
  )
}

# -----------------------------------------------------------------------------
# 3. DESeq2 and VST
# -----------------------------------------------------------------------------

dds_249 <- DESeq2::DESeqDataSetFromMatrix(
  countData = counts_249,
  colData = sample_info_249,
  design = ~ condition
)

dds_249 <- dds_249[
  rowSums(DESeq2::counts(dds_249)) >= 10,
]

dds_249 <- DESeq2::DESeq(dds_249)

vsd_249 <- DESeq2::vst(
  dds_249,
  blind = FALSE
)

saveRDS(
  dds_249,
  file.path(
    path_processed_249,
    "dds_GSE249311.rds"
  )
)

saveRDS(
  vsd_249,
  file.path(
    path_processed_249,
    "vst_GSE249311.rds"
  )
)

saveRDS(
  sample_info_249,
  file.path(
    path_processed_249,
    "sample_info_GSE249311.rds"
  )
)

# -----------------------------------------------------------------------------
# 4. Export prespecified disease-category contrasts
# -----------------------------------------------------------------------------

export_contrast <- function(
    dds_object,
    numerator,
    denominator = "Control"
) {
  result <- DESeq2::results(
    dds_object,
    contrast = c(
      "condition",
      numerator,
      denominator
    ),
    alpha = 0.05
  )

  result_df <- as.data.frame(result) |>
    tibble::rownames_to_column("gene") |>
    dplyr::arrange(padj)

  outfile <- file.path(
    path_tables_249,
    paste0(
      "DESeq2_",
      numerator,
      "_vs_",
      denominator,
      ".csv"
    )
  )

  write.csv(
    result_df,
    outfile,
    row.names = FALSE
  )

  invisible(result_df)
}

res_GDMA1_vs_Control <- export_contrast(
  dds_249,
  "GDMA1"
)

res_GDMA2_vs_Control <- export_contrast(
  dds_249,
  "GDMA2"
)

res_T2DM_vs_Control <- export_contrast(
  dds_249,
  "T2DM"
)

write_session_info("04_GSE249311_differential_expression.R")

message("\nPhase 4 complete.\n")
