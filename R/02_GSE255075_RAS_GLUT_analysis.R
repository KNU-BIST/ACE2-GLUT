# =============================================================================
# Phase 2 — Targeted RAS–GLUT analysis of GSE255075
# Script: 02_GSE255075_RAS_GLUT_analysis.R
# =============================================================================
#
# Project
# -------
# Placental ACE2 and glucose transporter expression in gestational diabetes
# mellitus
#
# Purpose
# -------
# Perform hypothesis-driven analysis of placental GLUT/SLC2 and
# renin-angiotensin-system genes in the GSE255075 VST expression matrix.
#
# This phase generates gene-of-interest expression plots, exploratory
# Normal-vs-GDM statistics, structured RAS/GLUT heatmaps, and Pearson
# correlation matrices.
#
# Reproducibility boundary
# ------------------------
# GSE255075 contains only three Normal and three GDM samples. Gene-level
# secondary tests and correlation analyses are therefore exploratory and should
# not be interpreted as independent mechanistic validation.
#
# The nominal correlation symbols reproduce the historical visualization.
# BH-adjusted correlation P values are additionally exported for transparency.
#
# Inputs
# ------
# data/processed/GSE255075/vst_GSE255075.rds
# data/processed/GSE255075/sample_info_GSE255075.rds
#
# Outputs
# -------
# outputs/tables/GSE255075/GOI_expression_long.csv
# outputs/tables/GSE255075/GOI_Welch_ttest_results.csv
# outputs/tables/GSE255075/GLUT_RAS_correlations.csv
# figures/exploratory/GSE255075/GOI_boxplots/
# figures/exploratory/GSE255075/GSE255075_GLUT_RAS_heatmap.*
# figures/exploratory/GSE255075/GSE255075_RAS_heatmap.*
# figures/exploratory/GSE255075/GSE255075_GLUT_RAS_correlation_tile.png
# =============================================================================

if (!requireNamespace("here", quietly = TRUE)) {
  stop("Install package 'here' before running this script.")
}
source(here::here("R", "00_project_setup.R"))

vsd_file <- file.path(path_processed_255, "vst_GSE255075.rds")
sample_file <- file.path(path_processed_255, "sample_info_GSE255075.rds")

assert_object_exists(vsd_file, "GSE255075 VST object")
assert_object_exists(sample_file, "GSE255075 sample metadata")

vsd_255 <- readRDS(vsd_file)
sample_info_255 <- readRDS(sample_file)
expr_mat_255 <- SummarizedExperiment::assay(vsd_255)

control_samples <- rownames(sample_info_255)[sample_info_255$condition == "Normal"]
gdm_samples <- rownames(sample_info_255)[sample_info_255$condition == "GDM"]

message("\n===== GSE255075: targeted RAS–GLUT analysis =====\n")

# -----------------------------------------------------------------------------
# 1. Helper for retrieving a gene from SYMBOL or ENSEMBL rownames
# -----------------------------------------------------------------------------

get_expr_for_gene <- function(expr_mat, gene_symbol) {
  rn <- rownames(expr_mat)
  rn_clean <- sub("\\..*$", "", rn)

  if (gene_symbol %in% rn) {
    return(
      list(
        gene = gene_symbol,
        expr = as.numeric(expr_mat[gene_symbol, ]),
        samples = colnames(expr_mat)
      )
    )
  }

  map <- tryCatch(
    clusterProfiler::bitr(
      gene_symbol,
      fromType = "SYMBOL",
      toType = "ENSEMBL",
      OrgDb = org.Hs.eg.db
    ),
    error = function(e) NULL
  )

  if (is.null(map) || nrow(map) == 0) {
    return(NULL)
  }

  idx <- which(rn_clean %in% unique(map$ENSEMBL))

  if (length(idx) == 0) {
    return(NULL)
  }

  if (length(idx) > 1) {
    means <- rowMeans(
      expr_mat[idx, , drop = FALSE],
      na.rm = TRUE
    )
    idx <- idx[which.max(means)]
  }

  list(
    gene = gene_symbol,
    expr = as.numeric(expr_mat[idx, ]),
    samples = colnames(expr_mat)
  )
}

plot_gene_box <- function(
    gene_symbol,
    expr_vec,
    sample_names,
    sample_info
) {
  df <- data.frame(
    sample = sample_names,
    expression = expr_vec
  ) |>
    dplyr::left_join(
      sample_info |>
        as.data.frame() |>
        tibble::rownames_to_column("sample"),
      by = "sample"
    )

  ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = condition,
      y = expression
    )
  ) +
    ggplot2::geom_boxplot(
      fill = "white",
      color = "black",
      linewidth = box_line,
      width = 0.6,
      outlier.shape = NA
    ) +
    ggplot2::geom_point(
      size = 2,
      position = ggplot2::position_jitter(
        width = 0.08,
        height = 0
      )
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::labs(
      title = gene_symbol,
      x = NULL,
      y = "VST-normalized expression"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        hjust = 0.5
      )
    )
}

# -----------------------------------------------------------------------------
# 2. GOI expression plots and exploratory Welch tests
# -----------------------------------------------------------------------------

goi_plot_dir <- file.path(path_fig_255, "GOI_boxplots")
dir.create(goi_plot_dir, recursive = TRUE, showWarnings = FALSE)

goi_long_255 <- dplyr::bind_rows(
  lapply(
    genes_of_interest,
    function(g) {
      hit <- get_expr_for_gene(expr_mat_255, g)

      if (is.null(hit)) {
        message("Skipping ", g, ": not found.")
        return(NULL)
      }

      p <- plot_gene_box(
        hit$gene,
        hit$expr,
        hit$samples,
        sample_info_255
      )

      save_plot(
        p,
        file.path(
          goi_plot_dir,
          paste0(safe_name(g), "_VST_boxplot.png")
        ),
        width = 4.6,
        height = 4.2
      )

      data.frame(
        sample = hit$samples,
        gene = g,
        expression = hit$expr
      )
    }
  )
) |>
  dplyr::left_join(
    sample_info_255 |>
      as.data.frame() |>
      tibble::rownames_to_column("sample"),
    by = "sample"
  ) |>
  dplyr::filter(!is.na(condition))

write.csv(
  goi_long_255,
  file.path(path_tables_255, "GOI_expression_long.csv"),
  row.names = FALSE
)

goi_ttests_255 <- goi_long_255 |>
  dplyr::group_by(gene) |>
  dplyr::group_modify(
    ~ {
      normal_values <- .x$expression[.x$condition == "Normal"]
      gdm_values <- .x$expression[.x$condition == "GDM"]

      tt <- tryCatch(
        stats::t.test(
          gdm_values,
          normal_values,
          var.equal = FALSE
        ),
        error = function(e) NULL
      )

      data.frame(
        n_Normal = length(normal_values),
        n_GDM = length(gdm_values),
        mean_Normal = mean(normal_values, na.rm = TRUE),
        mean_GDM = mean(gdm_values, na.rm = TRUE),
        diff_GDM_minus_Normal =
          mean(gdm_values, na.rm = TRUE) -
          mean(normal_values, na.rm = TRUE),
        p_raw = if (is.null(tt)) NA_real_ else tt$p.value
      )
    }
  ) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    p_BH = p.adjust(p_raw, method = "BH")
  ) |>
  dplyr::arrange(p_BH)

write.csv(
  goi_ttests_255,
  file.path(path_tables_255, "GOI_Welch_ttest_results.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 3. Structured GLUT + RAS heatmap
# -----------------------------------------------------------------------------

genes_present_255 <- intersect(
  gene_panel,
  rownames(expr_mat_255)
)

mat_255 <- as.matrix(
  expr_mat_255[
    genes_present_255,
    c(control_samples, gdm_samples),
    drop = FALSE
  ]
)

mode(mat_255) <- "numeric"

mat_scaled_255 <- t(scale(t(mat_255)))
mat_scaled_255[is.na(mat_scaled_255)] <- 0

gene_order_255 <- intersect(
  c(
    glut_genes,
    "REN", "AGT", "ACE", "ACE2", "AGTR1",
    "MME", "THOP1", "ANPEP", "LNPEP", "ADAM17", "TIMP3"
  ),
  rownames(mat_scaled_255)
)

mat_scaled_255 <- mat_scaled_255[
  gene_order_255,
  ,
  drop = FALSE
]

gene_group_255 <- ifelse(
  rownames(mat_scaled_255) %in% glut_genes,
  "GLUT transporters",
  ifelse(
    rownames(mat_scaled_255) %in% ras_core,
    "RAS core",
    "RAS regulators"
  )
)

cond_255 <- c(
  rep("Normal", length(control_samples)),
  rep("GDM", length(gdm_samples))
)

ha_255 <- ComplexHeatmap::HeatmapAnnotation(
  Condition = cond_255,
  col = list(
    Condition = c(
      Normal = "#4DBBD5",
      GDM = "#E64B35"
    )
  )
)

col_fun_255 <- circlize::colorRamp2(
  c(-2, 0, 2),
  c("#3B4CC0", "white", "#B40426")
)

ht_glut_ras_255 <- ComplexHeatmap::Heatmap(
  mat_scaled_255,
  name = "Z-score",
  col = col_fun_255,
  top_annotation = ha_255,
  row_split = gene_group_255,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_title = "GLUT and RAS gene expression — GSE255075",
  column_names_rot = 45
)

pdf(
  file.path(path_fig_255, "GSE255075_GLUT_RAS_heatmap.pdf"),
  width = 6,
  height = 7
)
ComplexHeatmap::draw(ht_glut_ras_255)
dev.off()

png(
  file.path(path_fig_255, "GSE255075_GLUT_RAS_heatmap.png"),
  width = 6,
  height = 7,
  units = "in",
  res = dpi_out
)
ComplexHeatmap::draw(ht_glut_ras_255)
dev.off()

# -----------------------------------------------------------------------------
# 4. RAS-only heatmap
# -----------------------------------------------------------------------------

ras_genes_present_255 <- intersect(
  ras_genes,
  rownames(expr_mat_255)
)

ras_mat_255 <- as.matrix(
  expr_mat_255[
    ras_genes_present_255,
    c(control_samples, gdm_samples),
    drop = FALSE
  ]
)

mode(ras_mat_255) <- "numeric"

ras_scaled_255 <- t(scale(t(ras_mat_255)))
ras_scaled_255[is.na(ras_scaled_255)] <- 0

ras_groups_255 <- ifelse(
  rownames(ras_scaled_255) %in% ras_core,
  "RAS core",
  "RAS regulators"
)

ht_ras_255 <- ComplexHeatmap::Heatmap(
  ras_scaled_255,
  name = "Z-score",
  col = col_fun_255,
  top_annotation = ha_255,
  row_split = ras_groups_255,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_title = "RAS pathway gene expression — GSE255075",
  column_names_rot = 45
)

pdf(
  file.path(path_fig_255, "GSE255075_RAS_heatmap.pdf"),
  width = 6,
  height = 7
)
ComplexHeatmap::draw(ht_ras_255)
dev.off()

png(
  file.path(path_fig_255, "GSE255075_RAS_heatmap.png"),
  width = 6,
  height = 7,
  units = "in",
  res = dpi_out
)
ComplexHeatmap::draw(ht_ras_255)
dev.off()

# -----------------------------------------------------------------------------
# 5. GLUT–RAS Pearson correlations
# -----------------------------------------------------------------------------

present_cor_255 <- intersect(
  c(glut_genes, reg_genes_cor),
  rownames(expr_mat_255)
)

expr_sub_cor_255 <- expr_mat_255[
  present_cor_255,
  ,
  drop = FALSE
]

expr_glut_255 <- t(
  expr_sub_cor_255[
    intersect(glut_genes, rownames(expr_sub_cor_255)),
    ,
    drop = FALSE
  ]
)

expr_reg_255 <- t(
  expr_sub_cor_255[
    intersect(reg_genes_cor, rownames(expr_sub_cor_255)),
    ,
    drop = FALSE
  ]
)

pairs_cor_255 <- expand.grid(
  GLUT = colnames(expr_glut_255),
  Reg = colnames(expr_reg_255),
  stringsAsFactors = FALSE
)

cor_list_255 <- lapply(
  seq_len(nrow(pairs_cor_255)),
  function(i) {
    ct <- stats::cor.test(
      expr_glut_255[, pairs_cor_255$GLUT[i]],
      expr_reg_255[, pairs_cor_255$Reg[i]],
      method = "pearson"
    )

    data.frame(
      GLUT = pairs_cor_255$GLUT[i],
      Reg = pairs_cor_255$Reg[i],
      r = unname(ct$estimate),
      p_raw = ct$p.value,
      stringsAsFactors = FALSE
    )
  }
)

cor_df_255 <- dplyr::bind_rows(cor_list_255) |>
  dplyr::mutate(
    p_BH = p.adjust(p_raw, method = "BH"),
    sign_nominal = cut(
      p_raw,
      breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
      labels = c("***", "**", "*", "")
    ),
    label = paste0(
      sprintf("%.2f", r),
      sign_nominal
    )
  )

write.csv(
  cor_df_255,
  file.path(path_tables_255, "GLUT_RAS_correlations.csv"),
  row.names = FALSE
)

plot_cor_255 <- cor_df_255 |>
  dplyr::mutate(
    GLUT = factor(GLUT, levels = rev(glut_genes)),
    Reg = factor(Reg, levels = reg_genes_cor)
  )

p_cor_255 <- ggplot2::ggplot(
  plot_cor_255,
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
    "GLUT–RAS correlations — GSE255075"
  )

save_plot(
  p_cor_255,
  file.path(
    path_fig_255,
    "GSE255075_GLUT_RAS_correlation_tile.png"
  ),
  width = 25 / 2.54,
  height = 18 / 2.54
)

write_session_info("02_GSE255075_RAS_GLUT_analysis.R")

message("\nPhase 2 complete.\n")
