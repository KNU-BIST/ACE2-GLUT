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
# Author
# ------
# Md. Masudul Haque
#
# Purpose
# -------
# Reproduce the historical GSE249311 placental RNA-seq analysis while using
# clinically accurate source-dataset terminology.
#
# Source dataset categories
# -------------------------
# Control = 11
# GDMA1   = 5   lifestyle-managed gestational diabetes
# GDMA2   = 9   medication-requiring gestational diabetes
# T2DM    = 5   pregestational type 2 diabetes
#
# IMPORTANT:
# GDMA2 must not be interpreted as synonymous with insulin-treated GDM.
# Aspirin exposure is unavailable in this public dataset.
#
# Historical reproducibility choices
# ----------------------------------
# Duplicate source feature labels are retained as separate rows with
# make.unique(), matching the historical analysis.
#
# Low-count filtering:
#   retain genes with total count > 10
#
# VST:
#   blind = TRUE, matching the historical downstream-expression analysis
#
# DESeq2:
#   design = ~ condition
#   Benjamini-Hochberg multiple-testing correction
#
# Historical primary comparison:
#   GDMA2 versus Control
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
#
# outputs/tables/GSE249311/
#   DESeq2_GDMA1_vs_Control.csv
#   DESeq2_GDMA2_vs_Control.csv
#   DESeq2_T2DM_vs_Control.csv
#   GSE249311_feature_id_map.csv
# =============================================================================


# -----------------------------------------------------------------------------
# 0. Project setup
# -----------------------------------------------------------------------------

if (!requireNamespace("here", quietly = TRUE)) {
  stop("Install package 'here' before running this script.")
}

source(
  here::here(
    "R",
    "00_project_setup.R"
  )
)


# -----------------------------------------------------------------------------
# 1. Input
# -----------------------------------------------------------------------------

input_xlsx_249 <- file.path(
  path_data_raw_249,
  "GSE249311_count_matrix.xlsx"
)

assert_file_exists(
  input_xlsx_249
)

message(
  "\n===== GSE249311: differential expression =====\n"
)

counts_df_249 <- readxl::read_excel(
  input_xlsx_249
)

if (ncol(counts_df_249) < 2) {
  stop(
    "GSE249311 count matrix must contain one gene column plus sample columns."
  )
}


# -----------------------------------------------------------------------------
# 2. Gene identifiers
# -----------------------------------------------------------------------------

gene_ids_249 <- trimws(
  as.character(
    counts_df_249[[1]]
  )
)

valid_gene_249 <- !is.na(gene_ids_249) &
  gene_ids_249 != ""

counts_df_249 <- counts_df_249[
  valid_gene_249,
  ,
  drop = FALSE
]

gene_ids_249 <- gene_ids_249[
  valid_gene_249
]


# Historical analysis retained duplicate feature rows separately.
# make.unique() creates deterministic suffixes such as ".1", ".2", etc.

analysis_gene_ids_249 <- make.unique(
  gene_ids_249
)

feature_map_249 <- data.frame(
  source_gene_id = gene_ids_249,
  analysis_gene_id = analysis_gene_ids_249,
  stringsAsFactors = FALSE
)

write.csv(
  feature_map_249,
  file.path(
    path_tables_249,
    "GSE249311_feature_id_map.csv"
  ),
  row.names = FALSE
)

n_duplicate_rows_249 <- sum(
  duplicated(gene_ids_249)
)

n_duplicate_ids_249 <- length(
  unique(
    gene_ids_249[
      duplicated(gene_ids_249) |
        duplicated(
          gene_ids_249,
          fromLast = TRUE
        )
    ]
  )
)

message(
  "Duplicate feature identifiers retained using make.unique(): ",
  n_duplicate_rows_249,
  " duplicate occurrences beyond the first across ",
  n_duplicate_ids_249,
  " source identifiers."
)


# -----------------------------------------------------------------------------
# 3. Count matrix
# -----------------------------------------------------------------------------

count_data_249 <- counts_df_249[
  ,
  -1,
  drop = FALSE
]

count_data_249[] <- lapply(
  count_data_249,
  function(x) {
    suppressWarnings(
      as.numeric(
        as.character(x)
      )
    )
  }
)

if (anyNA(count_data_249)) {

  bad_columns_249 <- names(count_data_249)[
    vapply(
      count_data_249,
      anyNA,
      logical(1)
    )
  ]

  stop(
    "Missing/non-numeric count values found in: ",
    paste(
      bad_columns_249,
      collapse = ", "
    )
  )
}

counts_249 <- as.matrix(
  count_data_249
)

rownames(counts_249) <- analysis_gene_ids_249

if (any(counts_249 < 0)) {
  stop(
    "Negative count values found."
  )
}

if (
  any(
    abs(
      counts_249 -
      round(counts_249)
    ) > 1e-8
  )
) {
  stop(
    "Non-integer count values detected in GSE249311."
  )
}

storage.mode(counts_249) <- "integer"

sample_names_249 <- colnames(
  counts_249
)

if (length(sample_names_249) != 30) {
  stop(
    "Expected 30 GSE249311 samples but found ",
    length(sample_names_249),
    "."
  )
}


# -----------------------------------------------------------------------------
# 4. Sample metadata
# -----------------------------------------------------------------------------
#
# Derive clinical categories from the actual sample names rather than relying
# only on column position. This is safer while reproducing the same grouping.

condition_249 <- ifelse(
  grepl("^Control_", sample_names_249),
  "Control",
  ifelse(
    grepl("^A1_", sample_names_249),
    "GDMA1",
    ifelse(
      grepl("^A2_", sample_names_249),
      "GDMA2",
      ifelse(
        grepl("^T2DM_", sample_names_249),
        "T2DM",
        NA_character_
      )
    )
  )
)

if (anyNA(condition_249)) {
  stop(
    "Unable to assign one or more GSE249311 samples to a clinical category."
  )
}

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

expected_counts_249 <- c(
  Control = 11,
  GDMA1 = 5,
  GDMA2 = 9,
  T2DM = 5
)

observed_counts_249 <- table(
  sample_info_249$condition
)

print(
  observed_counts_249
)

if (
  !identical(
    as.integer(observed_counts_249),
    as.integer(expected_counts_249)
  )
) {
  stop(
    "GSE249311 group counts do not match the expected 11/5/9/5 structure."
  )
}

write.csv(
  sample_info_249,
  file.path(
    path_tables_249,
    "sample_info_GSE249311.csv"
  ),
  row.names = TRUE
)


# -----------------------------------------------------------------------------
# 5. Construct DESeq2 dataset
# -----------------------------------------------------------------------------

dds_249 <- DESeq2::DESeqDataSetFromMatrix(
  countData = counts_249,
  colData = sample_info_249,
  design = ~ condition
)

# Historical filter: STRICTLY greater than 10.
dds_249 <- dds_249[
  rowSums(
    DESeq2::counts(dds_249)
  ) > 10,
]

message(
  "Genes retained after total-count >10 filter: ",
  nrow(dds_249)
)

if (nrow(dds_249) != 27536) {
  warning(
    "Expected 27,536 genes after historical >10 filter, but found ",
    nrow(dds_249),
    "."
  )
}


# -----------------------------------------------------------------------------
# 6. VST
# -----------------------------------------------------------------------------
#
# Historical downstream heatmaps, correlations and expression plots were based
# on the default blind VST. Set blind=TRUE explicitly for reproducibility.

vsd_249 <- DESeq2::vst(
  dds_249,
  blind = TRUE
)


# -----------------------------------------------------------------------------
# 7. DESeq2 model
# -----------------------------------------------------------------------------

dds_249 <- DESeq2::DESeq(
  dds_249
)


# -----------------------------------------------------------------------------
# 8. Save processed objects
# -----------------------------------------------------------------------------

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
# 9. Differential-expression contrasts
# -----------------------------------------------------------------------------

export_contrast_249 <- function(
    dds_object,
    numerator,
    denominator = "Control"
) {

  result_249 <- DESeq2::results(
    dds_object,
    contrast = c(
      "condition",
      numerator,
      denominator
    ),
    alpha = 0.1
  )

  result_df_249 <- as.data.frame(
    result_249
  ) |>
    tibble::rownames_to_column(
      "gene"
    ) |>
    dplyr::filter(
      !is.na(padj)
    ) |>
    dplyr::arrange(
      padj
    )

  outfile_249 <- file.path(
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
    result_df_249,
    outfile_249,
    row.names = FALSE
  )

  invisible(
    result_df_249
  )
}


# GDMA1: lifestyle-managed GDM versus Control
res_GDMA1_vs_Control <- export_contrast_249(
  dds_249,
  "GDMA1"
)


# GDMA2: medication-requiring GDM versus Control
# This is the terminology-corrected equivalent of the historical
# "GDM_insulin vs Control" contrast.
res_GDMA2_vs_Control <- export_contrast_249(
  dds_249,
  "GDMA2"
)


# Pregestational type 2 diabetes versus Control
res_T2DM_vs_Control <- export_contrast_249(
  dds_249,
  "T2DM"
)


# -----------------------------------------------------------------------------
# 10. Key reproducibility summary
# -----------------------------------------------------------------------------

message(
  "\nGDMA2 vs Control:"
)

message(
  "  FDR <0.05: ",
  sum(
    res_GDMA2_vs_Control$padj < 0.05,
    na.rm = TRUE
  )
)

message(
  "  FDR <0.05 and |log2FC| >1: ",
  sum(
    res_GDMA2_vs_Control$padj < 0.05 &
      abs(
        res_GDMA2_vs_Control$log2FoldChange
      ) > 1,
    na.rm = TRUE
  )
)


# -----------------------------------------------------------------------------
# 11. Session information
# -----------------------------------------------------------------------------

write_session_info(
  "04_GSE249311_differential_expression.R"
)

message(
  "\nPhase 4 complete.\n"
)