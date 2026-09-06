# =============================================================================
# Phase 5 — Targeted RAS–GLUT analysis of GSE249311
# Script: 05_GSE249311_RAS_GLUT_analysis.R
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
# Use the persisted GSE249311 DESeq2/VST objects to examine RAS–GLUT
# expression patterns, nominal Pearson correlations, GLUT expression
# distributions, DESeq2 Wald-test contrasts, ACE2 expression across clinical
# categories, and pathway context.
#
# Reproducibility boundary
# ------------------------
# GSE249311 is used as independent transcriptomic context.
#
# The source categories are Control, GDMA1, GDMA2, and T2DM. GDMA2 denotes
# medication-requiring GDM and is not interpreted as synonymous with insulin
# treatment. Aspirin exposure is unavailable.
#
# Correlation-heatmap significance symbols use nominal unadjusted Pearson
# P values to reproduce the manuscript analysis. BH-adjusted values are also
# exported and should be used when assessing multiplicity.
#
# Historical SLC2 group comparisons are reproduced using two-sided Welch
# t-tests on VST-normalized expression, with BH correction jointly across
# the 12 SLC2 comparisons.
#
# DESeq2 Wald-test contrasts are exported separately as transcriptome-level
# differential-expression inference.
#
# Inputs
# ------
# data/processed/GSE249311/dds_GSE249311.rds
# data/processed/GSE249311/vst_GSE249311.rds
# data/processed/GSE249311/sample_info_GSE249311.rds
#
# Outputs
# -------
# outputs/tables/GSE249311/
#   GLUT_RAS_correlations.csv
#   GLUT_VST_expression_long.csv
#   GLUT_VST_Welch_ttests.csv
#   GLUT_DESeq2_contrasts.csv
#   ACE2_DESeq2_contrasts.csv
#   ACE2_SLC2A9_correlation.csv
#   KEGG_GDMA2_vs_Control.csv
#
# figures/exploratory/GSE249311/
# =============================================================================

if (!requireNamespace("here", quietly = TRUE)) {
  stop("Install package 'here' before running this script.")
}
source(here::here("R", "00_project_setup.R"))

dds_file <- file.path(
  path_processed_249,
  "dds_GSE249311.rds"
)

vsd_file <- file.path(
  path_processed_249,
  "vst_GSE249311.rds"
)

sample_file <- file.path(
  path_processed_249,
  "sample_info_GSE249311.rds"
)

assert_object_exists(
  dds_file,
  "GSE249311 DESeq2 object"
)

assert_object_exists(
  vsd_file,
  "GSE249311 VST object"
)

assert_object_exists(
  sample_file,
  "GSE249311 sample metadata"
)

dds_249 <- readRDS(dds_file)
vsd_249 <- readRDS(vsd_file)
sample_info_249 <- readRDS(sample_file)

expr_mat_249 <- SummarizedExperiment::assay(
  vsd_249
)

message("\n===== GSE249311: targeted RAS–GLUT analysis =====\n")

# -----------------------------------------------------------------------------
# 1. Structured GLUT + RAS heatmap
# -----------------------------------------------------------------------------

present_249 <- intersect(
  gene_panel,
  rownames(expr_mat_249)
)

expr_panel_249 <- expr_mat_249[
  present_249,
  ,
  drop = FALSE
]

expr_panel_scaled_249 <- t(
  scale(t(expr_panel_249))
)

expr_panel_scaled_249 <- expr_panel_scaled_249[
  complete.cases(expr_panel_scaled_249),
  ,
  drop = FALSE
]

row_order_249 <- c(
  intersect(
    glut_genes,
    rownames(expr_panel_scaled_249)
  ),
  intersect(
    ras_core,
    rownames(expr_panel_scaled_249)
  ),
  intersect(
    ras_reg,
    rownames(expr_panel_scaled_249)
  )
)

expr_panel_scaled_249 <- expr_panel_scaled_249[
  row_order_249,
  ,
  drop = FALSE
]

gene_cat_249 <- data.frame(
  Category = ifelse(
    row_order_249 %in% glut_genes,
    "GLUT",
    ifelse(
      row_order_249 %in% ras_core,
      "RAS_core",
      "RAS_regulator"
    )
  ),
  row.names = row_order_249
)

ann_col_249 <- data.frame(
  Condition = sample_info_249$condition,
  row.names = rownames(sample_info_249)
)

ann_colors_249 <- list(
  Category = c(
    GLUT = "#1b9e77",
    RAS_core = "#d95f02",
    RAS_regulator = "#7570b3"
  ),
  Condition = c(
    Control = "#4DAF4A",
    GDMA1 = "#E41A1C",
    GDMA2 = "#377EB8",
    T2DM = "#984EA3"
  )
)

gaps_row_249 <- c(
  sum(row_order_249 %in% glut_genes),
  sum(row_order_249 %in% c(glut_genes, ras_core))
)

gaps_row_249 <- gaps_row_249[
  gaps_row_249 > 0 &
    gaps_row_249 < nrow(expr_panel_scaled_249)
]

png(
  file.path(
    path_fig_249,
    "GSE249311_GLUT_RAS_structured_heatmap.png"
  ),
  width = 8,
  height = 7,
  units = "in",
  res = dpi_out
)

pheatmap::pheatmap(
  expr_panel_scaled_249,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_row = gene_cat_249,
  annotation_col = ann_col_249,
  annotation_colors = ann_colors_249,
  gaps_row = gaps_row_249,
  color = grDevices::colorRampPalette(
    c("navy", "white", "firebrick3")
  )(100),
  fontsize_row = 10,
  main = "Placental RAS–GLUT expression — GSE249311"
)

dev.off()

# -----------------------------------------------------------------------------
# 2. GLUT–RAS Pearson correlation matrix
# -----------------------------------------------------------------------------

present_cor_249 <- intersect(
  c(glut_genes, reg_genes_cor),
  rownames(expr_mat_249)
)

expr_sub_249 <- expr_mat_249[
  present_cor_249,
  ,
  drop = FALSE
]

expr_glut_249 <- t(
  expr_sub_249[
    intersect(glut_genes, rownames(expr_sub_249)),
    ,
    drop = FALSE
  ]
)

expr_reg_249 <- t(
  expr_sub_249[
    intersect(reg_genes_cor, rownames(expr_sub_249)),
    ,
    drop = FALSE
  ]
)

pairs_cor_249 <- expand.grid(
  GLUT = colnames(expr_glut_249),
  Reg = colnames(expr_reg_249),
  stringsAsFactors = FALSE
)

cor_list_249 <- lapply(
  seq_len(nrow(pairs_cor_249)),
  function(i) {
    ct <- stats::cor.test(
      expr_glut_249[, pairs_cor_249$GLUT[i]],
      expr_reg_249[, pairs_cor_249$Reg[i]],
      method = "pearson"
    )

    data.frame(
      GLUT = pairs_cor_249$GLUT[i],
      Reg = pairs_cor_249$Reg[i],
      r = unname(ct$estimate),
      p_raw = ct$p.value,
      stringsAsFactors = FALSE
    )
  }
)

cor_df_249 <- dplyr::bind_rows(
  cor_list_249
) |>
  dplyr::mutate(
    p_BH = p.adjust(
      p_raw,
      method = "BH"
    ),
    sign_nominal = cut(
      p_raw,
      breaks = c(
        -Inf,
        0.001,
        0.01,
        0.05,
        Inf
      ),
      labels = c(
        "***",
        "**",
        "*",
        ""
      )
    ),
    label = paste0(
      sprintf("%.2f", r),
      sign_nominal
    )
  )

write.csv(
  cor_df_249,
  file.path(
    path_tables_249,
    "GLUT_RAS_correlations.csv"
  ),
  row.names = FALSE
)

plot_cor_249 <- cor_df_249 |>
  dplyr::mutate(
    GLUT = factor(
      GLUT,
      levels = rev(glut_genes)
    ),
    Reg = factor(
      Reg,
      levels = reg_genes_cor
    )
  )

p_cor_249 <- ggplot2::ggplot(
  plot_cor_249,
  ggplot2::aes(
    x = Reg,
    y = GLUT,
    fill = r
  )
) +
  ggplot2::geom_tile(
    color = "white",
    linewidth = 0.5
  ) +
  ggplot2::geom_text(
    ggplot2::aes(label = label),
    size = 4
  ) +
  ggplot2::scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Pearson r"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1
    ),
    axis.title = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(
      hjust = 0.5,
      face = "bold"
    )
  ) +
  ggplot2::ggtitle(
    "GLUT–RAS correlations — GSE249311"
  )

save_plot(
  p_cor_249,
  file.path(
    path_fig_249,
    "GSE249311_GLUT_RAS_correlation_tile.png"
  ),
  width = 25 / 2.54,
  height = 18 / 2.54
)

# -----------------------------------------------------------------------------
# 3. GLUT VST boxplots
# -----------------------------------------------------------------------------

glut_present_249 <- intersect(
  glut_genes,
  rownames(expr_mat_249)
)

glut_wide_249 <- as.data.frame(
  t(
    expr_mat_249[
      glut_present_249,
      ,
      drop = FALSE
    ]
  )
) |>
  tibble::rownames_to_column("Sample")

metadata_249 <- data.frame(
  Sample = rownames(sample_info_249),
  Condition = as.character(
    sample_info_249$condition
  ),
  stringsAsFactors = FALSE
)

glut_wide_249 <- dplyr::left_join(
  glut_wide_249,
  metadata_249,
  by = "Sample"
)

if (anyNA(glut_wide_249$Condition)) {
  stop(
    "Some GSE249311 samples failed metadata matching."
  )
}

glut_wide_249$Condition <- factor(
  glut_wide_249$Condition,
  levels = c(
    "Control",
    "GDMA1",
    "GDMA2",
    "T2DM"
  )
)

glut_long_249 <- glut_wide_249 |>
  tidyr::pivot_longer(
    cols = dplyr::all_of(glut_present_249),
    names_to = "Gene",
    values_to = "VST_expression"
  )

write.csv(
  glut_long_249,
  file.path(
    path_tables_249,
    "GLUT_VST_expression_long.csv"
  ),
  row.names = FALSE
)

shape_values_249 <- c(
  Control = 16,
  GDMA1 = 17,
  GDMA2 = 15,
  T2DM = 18
)

plot_list_249 <- lapply(
  glut_present_249,
  function(gene_name) {
    gene_df <- glut_long_249 |>
      dplyr::filter(
        Gene == gene_name
      )

    p <- ggplot2::ggplot(
      gene_df,
      ggplot2::aes(
        x = Condition,
        y = VST_expression,
        shape = Condition
      )
    ) +
      ggplot2::geom_boxplot(
        fill = "white",
        color = "black",
        linewidth = box_line,
        outlier.shape = NA
      ) +
      ggplot2::geom_jitter(
        width = 0.15,
        height = 0,
        size = point_size
      ) +
      ggplot2::scale_shape_manual(
        values = shape_values_249,
        drop = FALSE
      ) +
      ggplot2::theme_bw() +
      ggplot2::labs(
        title = gene_name,
        y = "VST-normalized expression",
        x = NULL
      ) +
      ggplot2::theme(
        legend.position = "none",
        plot.title = ggplot2::element_text(
          hjust = 0.5,
          size = 16,
          face = "bold"
        ),
        axis.text.x = ggplot2::element_text(
          angle = 45,
          hjust = 1
        ),
        panel.grid = ggplot2::element_blank()
      )

    save_plot(
      p,
      file.path(
        path_fig_249,
        paste0(
          "GLUT_",
          gene_name,
          "_expression.png"
        )
      ),
      width = 5,
      height = 5
    )

    p
  }
)

names(plot_list_249) <- glut_present_249

# -----------------------------------------------------------------------------
# 4. Historical GLUT VST Welch tests
# -----------------------------------------------------------------------------
#
# These tests reproduce the historical exploratory analysis performed on
# VST-normalized expression values.
#
# Three comparisons are made for each of four SLC2 genes:
#   GDMA1 vs Control
#   GDMA2 vs Control
#   T2DM  vs Control
#
# BH correction is applied jointly across all 12 comparisons.
#
# These exploratory tests are retained for historical reproducibility.
# DESeq2 Wald-test contrasts are reported separately below.

glut_vst_welch_results_249 <- dplyr::bind_rows(
  lapply(
    glut_present_249,
    function(gene_name) {

      gene_df <- glut_long_249 |>
        dplyr::filter(
          Gene == gene_name
        )

      comparisons <- list(
        c("Control", "GDMA1"),
        c("Control", "GDMA2"),
        c("Control", "T2DM")
      )

      dplyr::bind_rows(
        lapply(
          comparisons,
          function(comp) {

            sub_df <- gene_df |>
              dplyr::filter(
                Condition %in% comp
              )

            tt <- stats::t.test(
              VST_expression ~ Condition,
              data = sub_df,
              var.equal = FALSE
            )

            means <- tapply(
              sub_df$VST_expression,
              sub_df$Condition,
              mean
            )

            data.frame(
              Gene = gene_name,
              Comparison = paste(
                comp[2],
                "vs",
                comp[1]
              ),
              n_Control = sum(
                sub_df$Condition == comp[1]
              ),
              n_Comparison = sum(
                sub_df$Condition == comp[2]
              ),
              mean_Control = unname(
                means[comp[1]]
              ),
              mean_Comparison = unname(
                means[comp[2]]
              ),
              difference_Comparison_minus_Control =
                unname(means[comp[2]] - means[comp[1]]),
              p_raw = tt$p.value,
              stringsAsFactors = FALSE
            )
          }
        )
      )
    }
  )
)

glut_vst_welch_results_249$p_BH <- p.adjust(
  glut_vst_welch_results_249$p_raw,
  method = "BH"
)

write.csv(
  glut_vst_welch_results_249,
  file.path(
    path_tables_249,
    "GLUT_VST_Welch_ttests.csv"
  ),
  row.names = FALSE
)

message(
  "Historical VST Welch tests with BH <0.05: ",
  sum(
    glut_vst_welch_results_249$p_BH < 0.05,
    na.rm = TRUE
  ),
  " / ",
  nrow(glut_vst_welch_results_249)
)

# -----------------------------------------------------------------------------
# 5. SLC2 DESeq2 Wald-test contrasts
# -----------------------------------------------------------------------------

deseq2_contrasts <- list(
  GDMA1_vs_Control = c(
    "condition",
    "GDMA1",
    "Control"
  ),
  GDMA2_vs_Control = c(
    "condition",
    "GDMA2",
    "Control"
  ),
  T2DM_vs_Control = c(
    "condition",
    "T2DM",
    "Control"
  )
)

extract_gene_contrast <- function(
    contrast_name,
    contrast_vector,
    genes
) {
  result <- DESeq2::results(
    dds_249,
    contrast = contrast_vector,
    alpha = 0.1
  )

  result_df <- as.data.frame(result) |>
    tibble::rownames_to_column("Gene") |>
    dplyr::filter(
      Gene %in% genes
    ) |>
    dplyr::mutate(
      Contrast = contrast_name
    )

  result_df
}

glut_deseq2_results <- dplyr::bind_rows(
  lapply(
    names(deseq2_contrasts),
    function(contrast_name) {
      extract_gene_contrast(
        contrast_name,
        deseq2_contrasts[[contrast_name]],
        glut_present_249
      )
    }
  )
) |>
  dplyr::select(
    Contrast,
    Gene,
    baseMean,
    log2FoldChange,
    lfcSE,
    stat,
    pvalue,
    padj
  )

write.csv(
  glut_deseq2_results,
  file.path(
    path_tables_249,
    "GLUT_DESeq2_contrasts.csv"
  ),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 5. ACE2 expression and DESeq2 contrasts
# -----------------------------------------------------------------------------

if ("ACE2" %in% rownames(expr_mat_249)) {
  ace2_df_249 <- data.frame(
    Sample = colnames(expr_mat_249),
    Condition = sample_info_249$condition,
    ACE2_VST = as.numeric(
      expr_mat_249["ACE2", ]
    )
  )

  p_ace2_249 <- ggplot2::ggplot(
    ace2_df_249,
    ggplot2::aes(
      x = Condition,
      y = ACE2_VST,
      fill = Condition
    )
  ) +
    ggplot2::geom_boxplot(
      alpha = 0.7,
      outlier.shape = NA
    ) +
    ggplot2::geom_jitter(
      width = 0.15,
      size = 2
    ) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      y = "ACE2 VST-normalized expression",
      x = NULL,
      title = "Placental ACE2 — GSE249311"
    ) +
    ggplot2::theme(
      legend.position = "none"
    )

  save_plot(
    p_ace2_249,
    file.path(
      path_fig_249,
      "GSE249311_ACE2_expression.png"
    ),
    width = 6,
    height = 5
  )

  ace2_deseq2_results <- dplyr::bind_rows(
    lapply(
      names(deseq2_contrasts),
      function(contrast_name) {
        extract_gene_contrast(
          contrast_name,
          deseq2_contrasts[[contrast_name]],
          "ACE2"
        )
      }
    )
  ) |>
    dplyr::select(
      Contrast,
      Gene,
      baseMean,
      log2FoldChange,
      lfcSE,
      stat,
      pvalue,
      padj
    )

  write.csv(
    ace2_deseq2_results,
    file.path(
      path_tables_249,
      "ACE2_DESeq2_contrasts.csv"
    ),
    row.names = FALSE
  )
}

# Explicit manuscript-relevant ACE2-SLC2A9 correlation

if (
  all(
    c("ACE2", "SLC2A9") %in%
    rownames(expr_mat_249)
  )
) {

  ace2_slc2a9_249 <- stats::cor.test(
    as.numeric(
      expr_mat_249["ACE2", ]
    ),
    as.numeric(
      expr_mat_249["SLC2A9", ]
    ),
    method = "pearson"
  )

  ace2_slc2a9_table_249 <- data.frame(
    Gene1 = "ACE2",
    Gene2 = "SLC2A9",
    n = length(
      expr_mat_249["ACE2", ]
    ),
    Pearson_r = unname(
      ace2_slc2a9_249$estimate
    ),
    p_two_sided = ace2_slc2a9_249$p.value,
    CI_lower = ace2_slc2a9_249$conf.int[1],
    CI_upper = ace2_slc2a9_249$conf.int[2]
  )

  write.csv(
    ace2_slc2a9_table_249,
    file.path(
      path_tables_249,
      "ACE2_SLC2A9_correlation.csv"
    ),
    row.names = FALSE
  )

  message(
    sprintf(
      paste0(
        "ACE2-SLC2A9: r = %.4f, ",
        "two-sided P = %.4f, n = %d"
      ),
      unname(
        ace2_slc2a9_249$estimate
      ),
      ace2_slc2a9_249$p.value,
      length(
        expr_mat_249["ACE2", ]
      )
    )
  )
}

# -----------------------------------------------------------------------------
# 6. KEGG context using GDMA2 vs Control DE result
# -----------------------------------------------------------------------------

de_file_249 <- file.path(
  path_tables_249,
  "DESeq2_GDMA2_vs_Control.csv"
)

if (file.exists(de_file_249)) {
  res_gdma2 <- read.csv(
    de_file_249,
    stringsAsFactors = FALSE
  )

  sig_genes_249 <- res_gdma2$gene[
    !is.na(res_gdma2$padj) &
      res_gdma2$padj < 0.05
  ]

  if (length(sig_genes_249) > 0) {
    gene_entrez_249 <- clusterProfiler::bitr(
      sig_genes_249,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    )

    if (nrow(gene_entrez_249) > 0) {
      kegg_res_249 <- clusterProfiler::enrichKEGG(
        gene = gene_entrez_249$ENTREZID,
        organism = "hsa"
      )

      write.csv(
        as.data.frame(kegg_res_249),
        file.path(
          path_tables_249,
          "KEGG_GDMA2_vs_Control.csv"
        ),
        row.names = FALSE
      )

      if (nrow(as.data.frame(kegg_res_249)) > 0) {
        p_kegg_249 <- enrichplot::dotplot(
          kegg_res_249,
          showCategory = 15
        ) +
          ggplot2::ggtitle(
            "KEGG enrichment — GDMA2 vs Control"
          )

        save_plot(
          p_kegg_249,
          file.path(
            path_fig_249,
            "GSE249311_KEGG_GDMA2_vs_Control.png"
          ),
          width = 8,
          height = 6
        )
      }
    }
  }
}

write_session_info("05_GSE249311_RAS_GLUT_analysis.R")

message("\nPhase 5 complete.\n")
