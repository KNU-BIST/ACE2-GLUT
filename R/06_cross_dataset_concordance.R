# =============================================================================
# Phase 6 — Cross-dataset RAS–GLUT concordance
# Script: 06_cross_dataset_concordance.R
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
# Compare GLUT–RAS Pearson-correlation estimates from GSE255075 and GSE249311
# using persisted result tables rather than in-memory objects.
#
# Reproducibility boundary
# ------------------------
# This script summarizes consistency in correlation direction and nominal
# significance across two independent public datasets. It does not establish a
# common causal regulatory mechanism.
#
# GSE255075 contains only six placentas; cross-dataset concordance should
# therefore be interpreted descriptively.
#
# Inputs
# ------
# outputs/tables/GSE255075/GLUT_RAS_correlations.csv
# outputs/tables/GSE249311/GLUT_RAS_correlations.csv
#
# Outputs
# -------
# outputs/tables/GSE249311/cross_dataset_GLUT_RAS_concordance.csv
# =============================================================================

if (!requireNamespace("here", quietly = TRUE)) {
  stop("Install package 'here' before running this script.")
}
source(here::here("R", "00_project_setup.R"))

file_255 <- file.path(
  path_tables_255,
  "GLUT_RAS_correlations.csv"
)

file_249 <- file.path(
  path_tables_249,
  "GLUT_RAS_correlations.csv"
)

assert_object_exists(
  file_255,
  "GSE255075 correlation table"
)

assert_object_exists(
  file_249,
  "GSE249311 correlation table"
)

cor_255 <- read.csv(
  file_255,
  stringsAsFactors = FALSE
)

cor_249 <- read.csv(
  file_249,
  stringsAsFactors = FALSE
)

concordance <- merge(
  cor_255[
    ,
    c(
      "GLUT",
      "Reg",
      "r",
      "p_raw",
      "p_BH"
    )
  ],
  cor_249[
    ,
    c(
      "GLUT",
      "Reg",
      "r",
      "p_raw",
      "p_BH"
    )
  ],
  by = c(
    "GLUT",
    "Reg"
  ),
  suffixes = c(
    "_255",
    "_249"
  )
)

concordance <- concordance |>
  dplyr::mutate(
    same_direction =
      sign(r_255) == sign(r_249),
    nominal_significant_both =
      p_raw_255 < 0.05 &
      p_raw_249 < 0.05,
    BH_significant_both =
      p_BH_255 < 0.05 &
      p_BH_249 < 0.05,
    same_direction_and_nominal =
      same_direction &
      nominal_significant_both
  )

write.csv(
  concordance,
  file.path(
    path_tables_249,
    "cross_dataset_GLUT_RAS_concordance.csv"
  ),
  row.names = FALSE
)

summary_text <- c(
  paste(
    "Pairs compared:",
    nrow(concordance)
  ),
  paste(
    "Pairs with same correlation direction:",
    sum(concordance$same_direction, na.rm = TRUE)
  ),
  paste(
    "Pairs nominally significant in both datasets:",
    sum(
      concordance$nominal_significant_both,
      na.rm = TRUE
    )
  ),
  paste(
    "Pairs BH-significant in both datasets:",
    sum(
      concordance$BH_significant_both,
      na.rm = TRUE
    )
  )
)

writeLines(
  summary_text,
  file.path(
    path_tables_249,
    "cross_dataset_GLUT_RAS_concordance_summary.txt"
  )
)

write_session_info("06_cross_dataset_concordance.R")

message("\nPhase 6 complete.\n")
