############################################################
# Combined RNA-seq Analysis Pipeline
# Datasets : GSE255075 (Normal vs GDM, n = 3+3)
#            GSE249311 (Control / GDM_diet / GDM_insulin / T2DM)
# Author   : Masudul Haque
# Lab      : Neurosurgery Research Lab, KNUH, Daegu, South Korea
# Purpose  : DESeq2 | VST | Volcano | KEGG/GAGE | Pathview |
#            GO ORA + GSEA | GLUT/RAS heatmaps | Correlation tiles |
#            Expression boxplots | Pairwise statistics
############################################################

# ─── 0.  PACKAGES ────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
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
  library(RColorBrewer)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(reshape2)
  library(stringr)
  library(AnnotationDbi)
  library(openxlsx)
})

set.seed(123)


# ─── 1.  USER CONFIGURATION ──────────────────────────────────────────────────

# ── GSE255075 ──
input_csv_255    <- "~/R/GDM/GSE255075/GSE255075_2.csv"
gene_col_255     <- "GeneID"
control_samples  <- c("Normal1", "Normal2", "Normal3")
gdm_samples      <- c("GDM1",    "GDM2",    "GDM3")
outdir_255       <- "GSE255075_outputs"

# ── GSE249311 ──
input_xlsx_249   <- "GSE249311_count_matrix.xlsx"
outdir_249       <- "GSE249311_outputs"
# condition vector must match column order of count matrix
condition_249    <- c(rep("Control", 11), rep("GDM_diet", 5),
                      rep("GDM_insulin", 9), rep("T2DM", 5))

# ── Figure aesthetics ──
box_line   <- 0.8          # boxplot border linewidth
point_size <- 3            # jitter point size
dpi_out    <- 300          # output resolution

# ── Helper: save ggplot ──
save_plot <- function(p, filepath, width = 8, height = 5, dpi = dpi_out) {
  ggsave(filename = filepath, plot = p, width = width, height = height, dpi = dpi)
}

# ── Helper: robust saving for ggplot OR gglist (gseaplot2) ──
save_any_plot <- function(p, filepath, width_in = 7.5, height_in = 5.5, dpi = dpi_out) {
  ext <- tolower(tools::file_ext(filepath))
  if (ext == "pdf") {
    pdf(filepath, width = width_in, height = height_in, onefile = TRUE)
    print(p); dev.off()
  } else {
    png(filepath, width = width_in * dpi, height = height_in * dpi, res = dpi)
    print(p); dev.off()
  }
}

# ── Helper: safe filenames ──
safe_name <- function(x, maxlen = 80) {
  x <- str_replace_all(x, "[^A-Za-z0-9_\\-]+", "_")
  x <- str_replace_all(x, "_+", "_")
  x <- str_replace_all(x, "^_|_$", "")
  substr(x, 1, maxlen)
}

# Create output directories
for (d in c(outdir_255, outdir_249)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}


# ─── 2.  SHARED GENE PANELS ──────────────────────────────────────────────────

glut_genes <- c("SLC2A1", "SLC2A3", "SLC2A4", "SLC2A9")
ras_core   <- c("ACE", "ACE2", "REN", "AGT", "AGTR1", "AGTR2", "MAS1")
ras_reg    <- c("MME", "THOP1", "ANPEP", "LNPEP", "ADAM17", "TIMP3")
ras_genes  <- c(ras_core, ras_reg)

# full panel used for GLUT/RAS heatmaps and correlations
gene_panel <- c(glut_genes, ras_core, ras_reg)

# genes of interest for expression plots and volcano labels
genes_of_interest <- unique(c(glut_genes, ras_genes))


# ══════════════════════════════════════════════════════════════════════════════
#  PART A  ·  GSE255075   Normal vs GDM  (n = 3 + 3)
# ══════════════════════════════════════════════════════════════════════════════

message("\n===== PART A: GSE255075 =====\n")

# ── A-01  Load & clean count data ────────────────────────────────────────────

raw_data <- read.csv(input_csv_255, header = TRUE, sep = ",",
                     check.names = FALSE, stringsAsFactors = FALSE)

# Remove empty/NA column names (common GEO/Excel artefact)
if (any(colnames(raw_data) == "" | is.na(colnames(raw_data)))) {
  message("Removing empty column names …")
  raw_data <- raw_data[, !(colnames(raw_data) == "" | is.na(colnames(raw_data))),
                       drop = FALSE]
}

stopifnot(gene_col_255 %in% colnames(raw_data))
stopifnot(all(c(control_samples, gdm_samples) %in% colnames(raw_data)))

raw_data <- raw_data[, c(control_samples, gdm_samples, gene_col_255), drop = FALSE]
raw_data <- raw_data[!is.na(raw_data[[gene_col_255]]) & raw_data[[gene_col_255]] != "", ]

# Strip ENSEMBL version suffix if present
if (any(grepl("^ENSG", raw_data[[gene_col_255]])))
  raw_data[[gene_col_255]] <- sub("\\..*$", "", raw_data[[gene_col_255]])

raw_data <- raw_data |>
  distinct(.data[[gene_col_255]], .keep_all = TRUE) |>
  mutate(across(all_of(c(control_samples, gdm_samples)), as.numeric))

# Sanity checks
count_vals <- as.matrix(raw_data[, c(control_samples, gdm_samples)])
if (any(is.na(count_vals)))
  warning("NA values present in count columns.")
if (any(count_vals < 0, na.rm = TRUE))
  stop("Negative counts detected. DESeq2 requires non-negative raw counts.")
if (mean(abs(count_vals - round(count_vals)) > 1e-6, na.rm = TRUE) > 0.05)
  warning(">5% non-integer values: may be normalised data (not raw counts).")

count_df_255 <- raw_data |> column_to_rownames(var = gene_col_255)
count_mat_255 <- as.matrix(count_df_255)
storage.mode(count_mat_255) <- "integer"

write.csv(
  data.frame(GeneID = rownames(count_mat_255), count_mat_255),
  file = file.path(outdir_255, "A01_cleaned_count_matrix.csv"),
  row.names = FALSE
)


# ── A-02  Sample metadata ────────────────────────────────────────────────────

count_mat_255 <- count_mat_255[, c(control_samples, gdm_samples), drop = FALSE]

sample_info_255 <- data.frame(
  row.names = colnames(count_mat_255),
  condition = factor(c(rep("Normal", length(control_samples)),
                       rep("GDM",    length(gdm_samples))),
                     levels = c("Normal", "GDM"))
)

write.csv(sample_info_255, file.path(outdir_255, "A02_sample_info.csv"), row.names = TRUE)


# ── A-03  DESeq2 differential expression ────────────────────────────────────

dds_255 <- DESeqDataSetFromMatrix(
  countData = count_mat_255,
  colData   = sample_info_255,
  design    = ~ condition
)
dds_255 <- dds_255[rowSums(counts(dds_255)) >= 10, ]
dds_255 <- DESeq(dds_255)

res_255 <- results(dds_255, contrast = c("condition", "GDM", "Normal"))

# LFC shrinkage (apeglm)
coef_name_255 <- grep("^condition_.*GDM.*vs.*Normal$",
                      resultsNames(dds_255), value = TRUE)
res_shrunk_255 <- lfcShrink(dds_255, coef = coef_name_255, type = "apeglm")

res_df_255 <- as.data.frame(res_shrunk_255) |>
  rownames_to_column("gene") |>
  arrange(padj)

write.csv(res_df_255,
          file.path(outdir_255, "A03_DESeq2_shrunk_GDM_vs_Normal.csv"),
          row.names = FALSE)

deg_n      <- sum(!is.na(res_df_255$padj) & res_df_255$padj < 0.05)
deg_n_fc   <- sum(!is.na(res_df_255$padj) & res_df_255$padj < 0.05 &
                    abs(res_df_255$log2FoldChange) > 1)
cat("DEG summary (GSE255075)\n",
    "  padj<0.05 :", deg_n, "\n",
    "  padj<0.05 & |log2FC|>1 :", deg_n_fc, "\n",
    file = file.path(outdir_255, "A03_DEG_summary.txt"))


# ── A-04  QC: VST + PCA ──────────────────────────────────────────────────────

vsd_255   <- vst(dds_255, blind = FALSE)
expr_mat_255 <- assay(vsd_255)          # genes × samples

pca_df_255   <- plotPCA(vsd_255, intgroup = "condition", returnData = TRUE)
pct_var_255  <- round(100 * attr(pca_df_255, "percentVar"))

p_pca_255 <- ggplot(pca_df_255, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(title = "PCA – VST normalized counts (GSE255075)",
       x = paste0("PC1: ", pct_var_255[1], "%"),
       y = paste0("PC2: ", pct_var_255[2], "%"))

save_plot(p_pca_255, file.path(outdir_255, "A04_PCA_VST.png"), 6.5, 5)


# ── A-05  Volcano plots ──────────────────────────────────────────────────────

# Primary (apeglm-shrunk, labelled GOI)
png(file.path(outdir_255, "A05_volcano_GOI_labelled.png"),
    width = 1800, height = 1500, res = 250)
EnhancedVolcano(
  res_df_255,
  lab        = res_df_255$gene,
  x          = "log2FoldChange",
  y          = "padj",
  selectLab  = intersect(genes_of_interest, res_df_255$gene),
  pCutoff    = 0.05,
  FCcutoff   = 1,
  xlim       = c(-5, 5),
  ylim       = c(0, 10),
  pointSize  = 2.0,
  labSize    = 3.5,
  axisLabSize = 14,
  titleLabSize = 16,
  title      = "Control vs GDM (GSE255075)"
)
dev.off()

# Alternative full-label volcano (unshrunk res for full −log10 range)
png(file.path(outdir_255, "A05_volcano_full.png"),
    width = 2000, height = 2000, res = 300)
EnhancedVolcano(
  res_255,
  lab           = rownames(res_255),
  x             = "log2FoldChange",
  y             = "padj",
  border        = "full",
  borderWidth   = 1.5,
  borderColour  = "black",
  gridlines.major = FALSE,
  gridlines.minor = FALSE,
  xlim          = c(-5, 5),
  ylim          = c(0, 30),
  title         = "Control vs GDM (GSE255075)"
)
dev.off()


# ── A-06  GOI expression boxplots ────────────────────────────────────────────

outdir_goi_255 <- file.path(outdir_255, "A06_GOI_boxplots")
dir.create(outdir_goi_255, showWarnings = FALSE, recursive = TRUE)

# Helper: retrieve VST vector for a gene symbol (SYMBOL or ENSEMBL rownames)
get_expr_for_gene <- function(expr_mat, gene_symbol) {
  rn       <- rownames(expr_mat)
  rn_clean <- sub("\\..*$", "", rn)
  
  # Direct SYMBOL match
  if (gene_symbol %in% rn) {
    return(list(gene    = gene_symbol,
                expr    = as.numeric(expr_mat[gene_symbol, ]),
                samples = colnames(expr_mat)))
  }
  
  # SYMBOL -> ENSEMBL mapping
  map <- tryCatch(
    bitr(gene_symbol, fromType = "SYMBOL", toType = "ENSEMBL", OrgDb = org.Hs.eg.db),
    error = function(e) NULL
  )
  if (is.null(map) || nrow(map) == 0) return(NULL)
  
  idx <- which(rn_clean %in% unique(map$ENSEMBL))
  if (length(idx) == 0) return(NULL)
  if (length(idx) > 1) {
    means <- rowMeans(expr_mat[idx, , drop = FALSE], na.rm = TRUE)
    idx   <- idx[which.max(means)]
  }
  
  list(gene    = gene_symbol,
       expr    = as.numeric(expr_mat[idx, ]),
       samples = colnames(expr_mat))
}

# Plotting function (used for both dataset A and dataset B GOI plots)
plot_gene_box <- function(gene_symbol, expr_vec, sample_names, si) {
  df <- data.frame(sample = sample_names, expression = expr_vec) |>
    left_join(si |> rownames_to_column("sample"), by = "sample")
  
  ggplot(df, aes(x = condition, y = expression)) +
    geom_boxplot(fill = "white", color = "black", linewidth = box_line,
                 width = 0.6, outlier.shape = NA) +
    geom_point(size = 2, position = position_jitter(width = 0.08, height = 0)) +
    theme_bw(base_size = 13) +
    labs(title = gene_symbol, x = NULL, y = "VST-normalised expression") +
    theme(plot.title  = element_text(face = "bold", hjust = 0.5),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 11))
}

for (g in genes_of_interest) {
  hit <- get_expr_for_gene(expr_mat_255, g)
  if (is.null(hit)) { message("Skipping ", g, ": not found."); next }
  p <- plot_gene_box(hit$gene, hit$expr, hit$samples, sample_info_255)
  fname <- paste0(safe_name(g), "_VST_boxplot.png")
  ggsave(file.path(outdir_goi_255, fname), p, width = 4.6, height = 4.2, dpi = dpi_out)
}
message("A-06 done: GOI boxplots saved to ", outdir_goi_255)


# ── A-07  GOI t-test statistics ──────────────────────────────────────────────

outdir_stats_255 <- file.path(outdir_255, "A07_GOI_stats")
dir.create(outdir_stats_255, showWarnings = FALSE, recursive = TRUE)

# Ensure sample_info_255 has a 'sample' column for joining
si_255_long <- sample_info_255 |> rownames_to_column("sample")

goi_long_255 <- bind_rows(
  lapply(genes_of_interest, function(g) {
    hit <- get_expr_for_gene(expr_mat_255, g)
    if (is.null(hit)) return(NULL)
    data.frame(sample = hit$samples, gene = g, expression = hit$expr)
  })
) |>
  left_join(si_255_long, by = "sample") |>
  filter(!is.na(condition))

write.csv(goi_long_255, file.path(outdir_stats_255, "GOI_expression_long.csv"),
          row.names = FALSE)

ttest_res_255 <- goi_long_255 |>
  group_by(gene) |>
  summarise(
    n_Normal  = sum(condition == "Normal"),
    n_GDM     = sum(condition == "GDM"),
    mean_Normal = mean(expression[condition == "Normal"], na.rm = TRUE),
    mean_GDM    = mean(expression[condition == "GDM"],    na.rm = TRUE),
    diff_GDM_minus_Normal = mean_GDM - mean_Normal,
    t_test_p = tryCatch(
      t.test(expression ~ condition, data = cur_data(), var.equal = FALSE)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) |>
  mutate(FDR_BH = p.adjust(t_test_p, method = "BH")) |>
  arrange(FDR_BH)

write.csv(ttest_res_255, file.path(outdir_stats_255, "GOI_ttest_results.csv"),
          row.names = FALSE)
print(ttest_res_255)


# ── A-08  GLUT + RAS ComplexHeatmap ─────────────────────────────────────────

genes_present_255 <- intersect(gene_panel, rownames(expr_mat_255))
mat_255 <- as.matrix(expr_mat_255[genes_present_255, ])
mode(mat_255) <- "numeric"
mat_scaled_255 <- t(scale(t(mat_255)))
mat_scaled_255[is.na(mat_scaled_255)] <- 0

# Reorder to predefined biological order (keep only what's present)
gene_order_255 <- intersect(
  c(glut_genes, "REN", "AGT", "ACE", "ACE2", "AGTR1",
    "MME", "THOP1", "ANPEP", "LNPEP", "ADAM17", "TIMP3"),
  rownames(mat_scaled_255)
)
mat_scaled_255 <- mat_scaled_255[gene_order_255,
                                 c(control_samples, gdm_samples)]

gene_groups_255 <- c(
  rep("GLUT transporters", length(intersect(glut_genes, gene_order_255))),
  rep("RAS core",          length(intersect(c("REN","AGT","ACE","ACE2","AGTR1"), gene_order_255))),
  rep("RAS regulators",    length(intersect(c("MME","THOP1","ANPEP","LNPEP","ADAM17","TIMP3"), gene_order_255)))
)

cond_255 <- c(rep("Normal", length(control_samples)), rep("GDM", length(gdm_samples)))

ha_255 <- HeatmapAnnotation(
  Condition = cond_255,
  col = list(Condition = c(Normal = "#4DBBD5", GDM = "#E64B35"))
)
col_fun_255 <- colorRamp2(c(-2, 0, 2), c("#3B4CC0", "white", "#B40426"))

ht_glut_ras_255 <- Heatmap(
  mat_scaled_255,
  name             = "Z-score",
  col              = col_fun_255,
  top_annotation   = ha_255,
  row_split        = gene_groups_255,
  cluster_rows     = FALSE,
  cluster_columns  = FALSE,
  column_title     = "GLUT and RAS gene expression – Normal vs GDM (GSE255075)",
  column_names_rot = 45
)

pdf(file.path(outdir_255, "A08_GLUT_RAS_heatmap.pdf"), width = 6, height = 7)
draw(ht_glut_ras_255); dev.off()

png(file.path(outdir_255, "A08_GLUT_RAS_heatmap.png"),
    width = 2400, height = 2800, res = dpi_out)
draw(ht_glut_ras_255); dev.off()


# ── A-09  RAS-only ComplexHeatmap ────────────────────────────────────────────

ras_genes_present_255 <- intersect(ras_genes, rownames(expr_mat_255))
ras_mat_255 <- as.matrix(expr_mat_255[ras_genes_present_255,
                                      c(control_samples, gdm_samples)])
mode(ras_mat_255) <- "numeric"
ras_scaled_255 <- t(scale(t(ras_mat_255)))
ras_scaled_255[is.na(ras_scaled_255)] <- 0

ras_groups_255 <- ifelse(rownames(ras_scaled_255) %in% ras_core,
                         "RAS core", "RAS regulators")

ht_ras_255 <- Heatmap(
  ras_scaled_255,
  name             = "Z-score",
  col              = col_fun_255,
  top_annotation   = ha_255,
  row_split        = ras_groups_255,
  cluster_rows     = FALSE,
  cluster_columns  = FALSE,
  column_title     = "RAS pathway gene expression – Normal vs GDM (GSE255075)",
  column_names_rot = 45
)

pdf(file.path(outdir_255, "A09_RAS_heatmap.pdf"), width = 6, height = 7)
draw(ht_ras_255); dev.off()

png(file.path(outdir_255, "A09_RAS_heatmap.png"),
    width = 2400, height = 2800, res = dpi_out)
draw(ht_ras_255); dev.off()


# ── A-10  GLUT–RAS correlation tile (ggplot) ─────────────────────────────────

reg_genes_cor <- c("REN", "AGT", "ACE", "ACE2", "AGTR1", "MME", "THOP1",
                   "ANPEP", "LNPEP", "ADAM17", "TIMP3", "MAS1")

genes_needed_cor <- unique(c(glut_genes, reg_genes_cor))
present_cor_255  <- intersect(genes_needed_cor, rownames(expr_mat_255))

expr_sub_cor <- expr_mat_255[present_cor_255, ]
expr_glut_255 <- t(expr_sub_cor[intersect(glut_genes,    rownames(expr_sub_cor)), ])
expr_reg_255  <- t(expr_sub_cor[intersect(reg_genes_cor, rownames(expr_sub_cor)), ])

# Compute Pearson r and p for each GLUT × RAS pair
pairs_cor_255 <- expand.grid(
  GLUT = colnames(expr_glut_255),
  Reg  = colnames(expr_reg_255),
  stringsAsFactors = FALSE
)

cor_list_255 <- lapply(seq_len(nrow(pairs_cor_255)), function(i) {
  ct <- cor.test(expr_glut_255[, pairs_cor_255$GLUT[i]],
                 expr_reg_255[,  pairs_cor_255$Reg[i]],
                 method = "pearson")
  data.frame(GLUT = pairs_cor_255$GLUT[i], Reg = pairs_cor_255$Reg[i],
             r = unname(ct$estimate), p = ct$p.value, stringsAsFactors = FALSE)
})
cor_df_255 <- do.call(rbind, cor_list_255)

cor_df_255$sign  <- cut(cor_df_255$p, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                        labels = c("***", "**", "*", ""))
cor_df_255$label <- paste0(sprintf("%.2f", cor_df_255$r), cor_df_255$sign)
cor_df_255$GLUT  <- factor(cor_df_255$GLUT, levels = rev(glut_genes))
cor_df_255$Reg   <- factor(cor_df_255$Reg,  levels = reg_genes_cor)

p_cor_255 <- ggplot(cor_df_255, aes(x = Reg, y = GLUT, fill = r)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = label), size = 4) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0, limits = c(-1, 1), name = "Pearson r") +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(),
        axis.text.x  = element_text(angle = 45, hjust = 1),
        axis.title   = element_blank(),
        plot.title   = element_text(hjust = 0.5, face = "bold")) +
  ggtitle("Correlation: GLUT vs ACE2/RAS genes (GSE255075)")

save_plot(p_cor_255, file.path(outdir_255, "A10_GLUT_RAS_correlation_tile.png"),
          width = 25 / 2.54, height = 18 / 2.54)

write.csv(cor_df_255, file.path(outdir_255, "A10_GLUT_RAS_correlation_table.csv"),
          row.names = FALSE)


# ── A-11  GAGE KEGG pathway enrichment ───────────────────────────────────────

gene_ids_255 <- res_df_255$gene
is_ensembl_255 <- mean(grepl("^ENSG", gene_ids_255)) > 0.5

fc_df_255 <- res_df_255 |>
  dplyr::select(gene, log2FoldChange) |>
  filter(!is.na(log2FoldChange))

if (is_ensembl_255) {
  map_255 <- bitr(fc_df_255$gene, fromType = "ENSEMBL",
                  toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  fc_mapped_255 <- fc_df_255 |>
    inner_join(map_255, by = c("gene" = "ENSEMBL")) |>
    distinct(ENTREZID, .keep_all = TRUE)
} else {
  map_255 <- bitr(fc_df_255$gene, fromType = "SYMBOL",
                  toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  fc_mapped_255 <- fc_df_255 |>
    inner_join(map_255, by = c("gene" = "SYMBOL")) |>
    distinct(ENTREZID, .keep_all = TRUE)
}

foldchanges_255 <- setNames(
  fc_mapped_255$log2FoldChange,
  fc_mapped_255$ENTREZID
)
foldchanges_255 <- sort(foldchanges_255, decreasing = TRUE)

write.csv(fc_mapped_255, file.path(outdir_255, "A11_foldchange_ENTREZ_mapping.csv"),
          row.names = FALSE)

kegg.hsa_255    <- kegg.gsets(species = "hsa")
kegg.sets.hsa   <- kegg.hsa_255$kg.sets[kegg.hsa_255$sigmet.idx]
keggres_255     <- gage(foldchanges_255, gsets = kegg.sets.hsa, same.dir = TRUE)

if (!is.null(keggres_255$greater))
  write.csv(as.data.frame(keggres_255$greater) |> rownames_to_column("pathway") |> arrange(p.val),
            file.path(outdir_255, "A11_GAGE_KEGG_upregulated.csv"), row.names = FALSE)
if (!is.null(keggres_255$less))
  write.csv(as.data.frame(keggres_255$less) |> rownames_to_column("pathway") |> arrange(p.val),
            file.path(outdir_255, "A11_GAGE_KEGG_downregulated.csv"), row.names = FALSE)


# ── A-12  Pathway membership overlap + Pathview ──────────────────────────────

goi_map_255    <- bitr(genes_of_interest, fromType = "SYMBOL",
                       toType = "ENTREZID", OrgDb = org.Hs.eg.db)
goi_entrez_255 <- unique(goi_map_255$ENTREZID)
entrez2sym_255 <- setNames(goi_map_255$SYMBOL, goi_map_255$ENTREZID)

hits_255_all <- lapply(names(kegg.sets.hsa), function(pw) {
  ov <- intersect(goi_entrez_255, kegg.sets.hsa[[pw]])
  if (length(ov) > 0) ov else NULL
})
names(hits_255_all) <- names(kegg.sets.hsa)
hits_255 <- Filter(Negate(is.null), hits_255_all)

hit_tbl_255 <- data.frame(
  pathway         = names(hits_255),
  n_genes_present = sapply(hits_255, length),
  genes_present   = sapply(hits_255, function(x) paste(entrez2sym_255[x], collapse = ","))
) |> arrange(desc(n_genes_present), pathway)

write.csv(hit_tbl_255, file.path(outdir_255, "A12_KEGG_pathways_containing_GOI.csv"),
          row.names = FALSE)

up_paths_255   <- if (!is.null(keggres_255$greater)) rownames(keggres_255$greater) else character(0)
down_paths_255 <- if (!is.null(keggres_255$less))    rownames(keggres_255$less)    else character(0)

write.csv(hit_tbl_255 |> filter(pathway %in% up_paths_255),
          file.path(outdir_255, "A12_GOI_pathways_upregulated.csv"), row.names = FALSE)
write.csv(hit_tbl_255 |> filter(pathway %in% down_paths_255),
          file.path(outdir_255, "A12_GOI_pathways_downregulated.csv"), row.names = FALSE)

deg_tbl_255 <- res_df_255 |> filter(!is.na(padj) & padj < 0.05)
writeLines(
  c("GOI that are significant DEGs (padj<0.05):",
    paste(intersect(genes_of_interest, deg_tbl_255$gene), collapse = ", ")),
  file.path(outdir_255, "A12_GOI_DEG_overlap.txt")
)

# Pathview
pv_dir_255 <- file.path(outdir_255, "A12_pathview")
if (!dir.exists(pv_dir_255)) dir.create(pv_dir_255, recursive = TRUE)

up_hit_tbl_255   <- hit_tbl_255 |> filter(pathway %in% up_paths_255)
down_hit_tbl_255 <- hit_tbl_255 |> filter(pathway %in% down_paths_255)

path_ids_255 <- unique(c(
  substr(head(rownames(keggres_255$greater), 20), 1, 8),
  substr(up_hit_tbl_255$pathway,   1, 8),
  substr(down_hit_tbl_255$pathway, 1, 8)
))
path_ids_255 <- path_ids_255[nzchar(path_ids_255)]

oldwd <- getwd(); setwd(pv_dir_255)
invisible(sapply(path_ids_255, function(pid) {
  try(pathview(gene.data = foldchanges_255, pathway.id = pid, species = "hsa",
               kegg.native = TRUE, same.layer = FALSE,
               out.suffix = "GDM_vs_Normal"), silent = TRUE)
}))
setwd(oldwd)
writeLines(path_ids_255, file.path(outdir_255, "A12_pathview_plotted_pathway_ids.txt"))


# ── A-13  GO ORA + GSEA ──────────────────────────────────────────────────────

outdir_go_255 <- file.path(outdir_255, "A13_GO")
dir.create(outdir_go_255, showWarnings = FALSE, recursive = TRUE)

# Map gene IDs to ENTREZ for GO
if (is_ensembl_255) {
  res_go_filtered <- res_df_255 |>
    filter(!is.na(log2FoldChange)) |>
    mutate(gene = sub("\\..*$", "", gene))
  map_go_255 <- bitr(res_go_filtered$gene, fromType = "ENSEMBL",
                     toType = c("ENTREZID", "SYMBOL"), OrgDb = org.Hs.eg.db)
  res_mapped_go_255 <- res_go_filtered |>
    inner_join(map_go_255, by = c("gene" = "ENSEMBL")) |>
    distinct(ENTREZID, .keep_all = TRUE)
} else {
  res_go_filtered <- res_df_255 |>
    filter(!is.na(log2FoldChange))
  map_go_255 <- bitr(res_go_filtered$gene, fromType = "SYMBOL",
                     toType = c("ENTREZID", "SYMBOL"), OrgDb = org.Hs.eg.db)
  res_mapped_go_255 <- res_go_filtered |>
    inner_join(map_go_255, by = c("gene" = "SYMBOL")) |>
    distinct(ENTREZID, .keep_all = TRUE)
}
write.csv(res_mapped_go_255, file.path(outdir_go_255, "01_res_mapped_to_ENTREZ.csv"),
          row.names = FALSE)

# ORA
padj_go <- 0.05; lfc_go <- 1
deg_all_255  <- res_mapped_go_255 |> filter(!is.na(padj) & padj < padj_go)
deg_up_255   <- deg_all_255 |> filter(log2FoldChange >=  lfc_go)
deg_down_255 <- deg_all_255 |> filter(log2FoldChange <= -lfc_go)

run_enrichGO <- function(entrez_vec, ont) {
  enrichGO(gene = entrez_vec, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
           ont = ont, pAdjustMethod = "BH",
           pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
}

ego_bp   <- run_enrichGO(deg_all_255$ENTREZID,  "BP")
ego_mf   <- run_enrichGO(deg_all_255$ENTREZID,  "MF")
ego_cc   <- run_enrichGO(deg_all_255$ENTREZID,  "CC")
ego_up   <- run_enrichGO(deg_up_255$ENTREZID,   "BP")
ego_down <- run_enrichGO(deg_down_255$ENTREZID, "BP")

for (nm in c("ego_bp","ego_mf","ego_cc","ego_up","ego_down")) {
  write.csv(as.data.frame(get(nm)),
            file.path(outdir_go_255, paste0("02_ORA_", nm, ".csv")),
            row.names = FALSE)
}

# ORA dotplots
for (nm in c("ego_bp","ego_up","ego_down")) {
  obj <- get(nm)
  lbl <- switch(nm, ego_bp = "BP all sig", ego_up = "BP up", ego_down = "BP down")
  p_dot <- dotplot(obj, showCategory = 20) +
    ggtitle(paste("GO ORA –", lbl, "(GSE255075)")) + theme_bw()
  save_plot(p_dot, file.path(outdir_go_255, paste0("03_dotplot_ORA_", nm, ".png")),
            9, 6)
}

# Simplify BP ORA
ego_bp_simpl <- tryCatch(
  simplify(ego_bp, cutoff = 0.7, by = "p.adjust", select_fun = min),
  error = function(e) { message("simplify() unavailable, skipping."); ego_bp }
)
write.csv(as.data.frame(ego_bp_simpl),
          file.path(outdir_go_255, "04_ORA_GO_BP_simplified.csv"), row.names = FALSE)
save_plot(
  dotplot(ego_bp_simpl, showCategory = 20) +
    ggtitle("GO BP (ORA) – simplified (GSE255075)") + theme_bw(),
  file.path(outdir_go_255, "04_dotplot_ORA_BP_simplified.png"), 9, 6)

# GSEA
geneList_255 <- sort(setNames(res_mapped_go_255$log2FoldChange,
                              res_mapped_go_255$ENTREZID), decreasing = TRUE)
gsea_BP_255  <- gseGO(geneList = geneList_255, OrgDb = org.Hs.eg.db,
                      keyType = "ENTREZID", ont = "BP",
                      minGSSize = 10, maxGSSize = 500,
                      pvalueCutoff = 0.05, pAdjustMethod = "BH",
                      verbose = FALSE)

write.csv(as.data.frame(gsea_BP_255),
          file.path(outdir_go_255, "05_GSEA_GO_BP.csv"), row.names = FALSE)


# ── A-14  GSEA dotplots + enrichment curves ───────────────────────────────────

outdir_gsea_255 <- file.path(outdir_255, "A14_GSEA_plots")
dir.create(outdir_gsea_255, showWarnings = FALSE, recursive = TRUE)

# NES-split dotplot (up vs down)
gsea_res_255 <- as.data.frame(gsea_BP_255)
gsea_res_255$Direction <- ifelse(gsea_res_255$NES >= 0, "Upregulated", "Downregulated")
gsea_BP_255b <- gsea_BP_255; gsea_BP_255b@result <- gsea_res_255

p_gsea_split <- dotplot(gsea_BP_255b, showCategory = 20, x = "NES",
                        label_format = 45) +
  facet_grid(Direction ~ ., scales = "free_y") +
  theme_classic(base_size = 12) +
  theme(axis.text  = element_text(color = "black"),
        strip.text = element_text(color = "black"),
        plot.title = element_text(color = "black")) +
  labs(title = "GO BP (GSEA): Up vs Down in GDM (GSE255075)", x = "NES", y = NULL)

save_plot(p_gsea_split,
          file.path(outdir_gsea_255, "A14_GSEA_BP_up_vs_down.png"),
          17 / 2.54, 18 / 2.54)

# Save as high-quality PDF (preferred for manuscripts)
ggsave(file.path(outdir_gsea_255, "A14_GSEA_BP_up_vs_down.pdf"),
       p_gsea_split, width = 17 / 2.54, height = 18 / 2.54)

# Per-term enrichment plots (upregulated, padj < 0.05)
up_terms_255 <- gsea_res_255 |>
  filter(!is.na(NES), NES > 0, !is.na(p.adjust), p.adjust <= 0.05) |>
  arrange(p.adjust)

if (nrow(up_terms_255) > 0) {
  q_col_255 <- if ("qvalues" %in% colnames(up_terms_255)) "qvalues" else
    if ("qvalue"  %in% colnames(up_terms_255)) "qvalue" else NULL
  
  for (i in seq_len(nrow(up_terms_255))) {
    tid  <- up_terms_255$ID[i]
    tdsc <- up_terms_255$Description[i]
    nes_v <- round(up_terms_255$NES[i], 2)
    pad_v <- signif(up_terms_255$p.adjust[i], 3)
    q_v   <- if (!is.null(q_col_255)) signif(up_terms_255[[q_col_255]][i], 3) else NA
    
    ttl <- paste0(tdsc, "\nNES = ", nes_v, " | adj.P = ", pad_v,
                  if (!is.na(q_v)) paste0(" | q(FDR) = ", q_v) else "")
    
    # pvalue_table = FALSE avoids the Annotated/aplot object that routes
    # rendering to the RStudio graphics window instead of the png device
    # (causing 0 KB output). NES and p-values are embedded in ttl above.
    p_enr <- tryCatch(
      gseaplot2(gsea_BP_255, geneSetID = tid, title = ttl, pvalue_table = FALSE),
      error = function(e) { message("Plot error: ", tid); NULL }
    )
    if (is.null(p_enr)) next
    
    # IMPORTANT: sanitize BOTH description AND GO term ID.
    # GO term IDs contain a colon (e.g. GO:0007565). On Windows, colons are
    # illegal in filenames — writing to a path with a colon silently creates
    # a 0-byte extensionless file and hides the content in an NTFS Alternate
    # Data Stream, which is invisible to Explorer and list.files().
    safe_tid  <- gsub("[^A-Za-z0-9_-]", "_", tid)
    safe_desc <- gsub("[^A-Za-z0-9_-]", "_", substr(tdsc, 1, 60))
    fname     <- paste0(sprintf("%03d", i), "_", safe_desc, "_", safe_tid, ".png")
    fpath     <- file.path(outdir_gsea_255, fname)
    
    ggsave(fpath, plot = p_enr, width = 7.5, height = 5.5, dpi = dpi_out)
    message(sprintf("[%d/%d] %d bytes -> %s",
                    i, nrow(up_terms_255), file.info(fpath)$size, basename(fpath)))
  }
  message("A-14 done: ", nrow(up_terms_255), " enrichment plots saved.")
} else {
  message("A-14: No upregulated GSEA terms with padj<0.05.")
}

# Session info for GSE255075 part
sink(file.path(outdir_255, "A_sessionInfo.txt"))
sessionInfo(); sink()

message("\n----- Part A (GSE255075) complete -----\n")


# ══════════════════════════════════════════════════════════════════════════════
#  PART B  ·  GSE249311   Control / GDM_diet / GDM_insulin / T2DM
# ══════════════════════════════════════════════════════════════════════════════

message("\n===== PART B: GSE249311 =====\n")

# ── B-01  Load count matrix (Excel) ──────────────────────────────────────────

counts_df_249 <- read_excel(input_xlsx_249)
gene_ids_249  <- counts_df_249[[1]]
counts_249    <- as.matrix(counts_df_249[, -1])
rownames(counts_249) <- make.unique(as.character(gene_ids_249))
mode(counts_249) <- "numeric"
sample_names_249 <- colnames(counts_249)


# ── B-02  Sample metadata ────────────────────────────────────────────────────

stopifnot(length(condition_249) == length(sample_names_249))

sample_info_249 <- data.frame(
  row.names = sample_names_249,
  condition = factor(condition_249)
)
print(table(sample_info_249$condition))


# ── B-03  DESeq2 normalization (VST) + DE ───────────────────────────────────

dds_249 <- DESeqDataSetFromMatrix(
  countData = counts_249,
  colData   = sample_info_249,
  design    = ~ condition
)
dds_249 <- dds_249[rowSums(counts(dds_249)) > 10, ]

vsd_249    <- vst(dds_249)
expr_mat_249 <- assay(vsd_249)   # genes × samples

# Full DE (GDM_insulin vs Control as primary contrast; others can be added)
dds_249 <- DESeq(dds_249)

res_249 <- results(dds_249, contrast = c("condition", "GDM_insulin", "Control"))
res_df_249 <- as.data.frame(res_249) |>
  rownames_to_column("gene") |>
  filter(!is.na(padj)) |>
  arrange(padj)
write.csv(res_df_249, file.path(outdir_249, "B03_DEG_GDMinsulin_vs_Control.csv"),
          row.names = FALSE)


# ── B-04  Structured GLUT+RAS heatmap (pheatmap) ────────────────────────────

present_249    <- intersect(gene_panel, rownames(expr_mat_249))
expr_panel_249 <- expr_mat_249[present_249, ]

# Row-scale
expr_panel_scaled_249 <- t(scale(t(expr_panel_249)))
expr_panel_scaled_249 <- expr_panel_scaled_249[complete.cases(expr_panel_scaled_249), ]

# Gene category annotation (only for genes that survived scaling)
gene_cat_249 <- data.frame(
  Category = c(
    rep("GLUT",          sum(rownames(expr_panel_scaled_249) %in% glut_genes)),
    rep("RAS_core",      sum(rownames(expr_panel_scaled_249) %in% ras_core)),
    rep("RAS_regulator", sum(rownames(expr_panel_scaled_249) %in% ras_reg))
  ),
  row.names = rownames(expr_panel_scaled_249)
)

# Re-order rows by category (GLUT → RAS_core → RAS_regulator)
row_order_249 <- c(
  intersect(glut_genes, rownames(expr_panel_scaled_249)),
  intersect(ras_core,   rownames(expr_panel_scaled_249)),
  intersect(ras_reg,    rownames(expr_panel_scaled_249))
)
expr_panel_scaled_249 <- expr_panel_scaled_249[row_order_249, ]
gene_cat_249          <- gene_cat_249[row_order_249, , drop = FALSE]

ann_col_249 <- data.frame(Condition = sample_info_249$condition,
                          row.names = rownames(sample_info_249))
ann_colors_249 <- list(
  Category  = c(GLUT = "#1b9e77", RAS_core = "#d95f02", RAS_regulator = "#7570b3"),
  Condition = c(Control = "#4DAF4A", GDM_diet = "#E41A1C",
                GDM_insulin = "#377EB8", T2DM = "#984EA3")
)
gaps_row_249 <- c(
  sum(rownames(gene_cat_249) %in% glut_genes),
  sum(rownames(gene_cat_249) %in% c(glut_genes, ras_core))
)

png(file.path(outdir_249, "B04_GLUT_RAS_structured_heatmap.png"),
    width = 2400, height = 2200, res = dpi_out)
pheatmap(
  expr_panel_scaled_249,
  cluster_rows       = FALSE,
  cluster_cols       = FALSE,
  annotation_row     = gene_cat_249,
  annotation_col     = ann_col_249,
  annotation_colors  = ann_colors_249,
  gaps_row           = gaps_row_249,
  color              = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  fontsize_row       = 10,
  main               = "Placental RAS–GLUT metabolic axis (GSE249311)"
)
dev.off()


# ── B-05  GLUT–RAS Pearson correlation tile ──────────────────────────────────

present_cor_249  <- intersect(c(glut_genes, reg_genes_cor), rownames(expr_mat_249))
expr_sub_249     <- expr_mat_249[present_cor_249, ]
expr_glut_249    <- t(expr_sub_249[intersect(glut_genes,    rownames(expr_sub_249)), ])
expr_reg_249     <- t(expr_sub_249[intersect(reg_genes_cor, rownames(expr_sub_249)), ])

pairs_cor_249 <- expand.grid(
  GLUT = colnames(expr_glut_249),
  Reg  = colnames(expr_reg_249),
  stringsAsFactors = FALSE
)

cor_list_249 <- lapply(seq_len(nrow(pairs_cor_249)), function(i) {
  ct <- cor.test(expr_glut_249[, pairs_cor_249$GLUT[i]],
                 expr_reg_249[,  pairs_cor_249$Reg[i]],
                 method = "pearson")
  data.frame(GLUT = pairs_cor_249$GLUT[i], Reg = pairs_cor_249$Reg[i],
             r = unname(ct$estimate), p = ct$p.value, stringsAsFactors = FALSE)
})
cor_df_249 <- do.call(rbind, cor_list_249)
cor_df_249$sign  <- cut(cor_df_249$p, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                        labels = c("***", "**", "*", ""))
cor_df_249$label <- paste0(sprintf("%.2f", cor_df_249$r), cor_df_249$sign)
cor_df_249$GLUT  <- factor(cor_df_249$GLUT, levels = rev(glut_genes))
cor_df_249$Reg   <- factor(cor_df_249$Reg,  levels = reg_genes_cor)

p_cor_249 <- ggplot(cor_df_249, aes(x = Reg, y = GLUT, fill = r)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = label), size = 4) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0, limits = c(-1, 1), name = "Pearson r") +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(),
        axis.text.x  = element_text(angle = 45, hjust = 1),
        axis.title   = element_blank(),
        plot.title   = element_text(hjust = 0.5, face = "bold")) +
  ggtitle("Correlation: GLUT vs ACE2/RAS genes (GSE249311)")

save_plot(p_cor_249, file.path(outdir_249, "B05_GLUT_RAS_correlation_tile.png"),
          25 / 2.54, 18 / 2.54)
write.csv(cor_df_249, file.path(outdir_249, "B05_GLUT_RAS_correlation_table.csv"),
          row.names = FALSE)


# ── B-06  GLUT expression boxplots ───────────────────────────────────────────

outdir_goi_249 <- file.path(outdir_249, "B06_GLUT_boxplots")
dir.create(outdir_goi_249, showWarnings = FALSE, recursive = TRUE)

glut_present_249 <- intersect(glut_genes, rownames(expr_mat_249))

# Build long-format table for boxplots and statistics
glut_df_249 <- as.data.frame(t(expr_mat_249[glut_present_249, ]))
glut_df_249$Condition <- sample_info_249$condition
glut_long_249 <- melt(glut_df_249, id.vars = "Condition",
                      variable.name = "Gene", value.name = "VST_expression")

shape_values_249 <- c(Control = 16, GDM_diet = 17, GDM_insulin = 15, T2DM = 18)

for (g in glut_present_249) {
  gene_df <- glut_long_249[glut_long_249$Gene == g, ]
  
  p_glut <- ggplot(gene_df,
                   aes(x = Condition, y = VST_expression, shape = Condition)) +
    geom_boxplot(fill = "white", color = "black", linewidth = box_line,
                 outlier.shape = NA) +
    geom_jitter(width = 0.15, size = point_size) +
    scale_shape_manual(values = shape_values_249) +
    theme_bw() +
    labs(title = g, y = "VST normalised expression", x = NULL) +
    theme(legend.position  = "none",
          plot.title       = element_text(hjust = 0.5, size = 16, face = "bold"),
          axis.text.y      = element_text(size = 14),
          axis.title.y     = element_text(size = 15),
          axis.text.x      = element_text(angle = 45, hjust = 1, size = 14),
          panel.grid       = element_blank())
  
  ggsave(file.path(outdir_goi_249, paste0("GLUT_", g, "_expression.png")),
         p_glut, width = 5, height = 5, dpi = dpi_out)
}


# ── B-07  Pairwise GLUT statistics (t-test vs Control, FDR) ─────────────────

pairwise_res_249 <- glut_long_249 |>
  group_by(Gene) |>
  do({
    df <- .
    comparisons <- list(c("Control", "GDM_diet"),
                        c("Control", "GDM_insulin"),
                        c("Control", "T2DM"))
    do.call(rbind, lapply(comparisons, function(comp) {
      sub_df <- df[df$Condition %in% comp, ]
      test   <- t.test(VST_expression ~ Condition, data = sub_df)
      data.frame(comparison = paste(comp, collapse = "_vs_"), p_value = test$p.value)
    }))
  }) |>
  ungroup()

pairwise_res_249$FDR <- p.adjust(pairwise_res_249$p_value, method = "fdr")

final_table_249 <- pairwise_res_249 |>
  pivot_wider(names_from = comparison, values_from = c(p_value, FDR))

print(final_table_249)

write.xlsx(final_table_249, file.path(outdir_249, "B07_GLUT_pairwise_stats.xlsx"))
write.csv(pairwise_res_249, file.path(outdir_249, "B07_GLUT_pairwise_stats_long.csv"),
          row.names = FALSE)


# ── B-08  KEGG enrichment (clusterProfiler) ──────────────────────────────────

sig_genes_249 <- res_df_249$gene[res_df_249$padj < 0.05]

if (length(sig_genes_249) > 0) {
  gene_entrez_249 <- bitr(sig_genes_249, fromType = "SYMBOL",
                          toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  kegg_res_249 <- enrichKEGG(gene = gene_entrez_249$ENTREZID, organism = "hsa")
  
  p_kegg_249 <- dotplot(kegg_res_249, showCategory = 15) +
    ggtitle("KEGG pathway enrichment: GDM_insulin vs Control (GSE249311)")
  save_plot(p_kegg_249, file.path(outdir_249, "B08_KEGG_enrichment.png"), 8, 6)
  write.csv(as.data.frame(kegg_res_249),
            file.path(outdir_249, "B08_KEGG_enrichment_results.csv"), row.names = FALSE)
} else {
  message("B-08: No significant DEGs for KEGG enrichment.")
}

# ACE2 expression trajectory across the glycemic gradient 
ace2_df <- data.frame(
  sample    = colnames(expr_mat_249),
  Condition = sample_info_249$condition,
  ACE2      = expr_mat_249["ACE2", ]
)
p_ace2 <- ggplot(ace2_df, aes(Condition, ACE2, fill = Condition)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2) +
  theme_bw() +
  labs(y = "ACE2 VST expression", x = NULL,
       title = "Placental ACE2 across glycemic spectrum (GSE249311)")
ggsave(file.path(outdir_249, "B09_ACE2_trajectory.png"),
       p_ace2, width = 6, height = 5, dpi = 300)

# Pairwise tests vs Control
ace2_stats <- lapply(c("GDM_diet", "GDM_insulin", "T2DM"), function(g) {
  sub <- ace2_df[ace2_df$Condition %in% c("Control", g), ]
  t.test(ACE2 ~ Condition, data = sub) |>
    (\(t) data.frame(group = g, p = t$p.value, mean_diff = -diff(t$estimate)))()
})
ace2_stats_df <- do.call(rbind, ace2_stats)
ace2_stats_df$FDR <- p.adjust(ace2_stats_df$p, "BH")
write.csv(ace2_stats_df, file.path(outdir_249, "B09_ACE2_pairwise_stats.csv"),
          row.names = FALSE)



# full RAS panel (ACE, AGT, AGTR1, ADAM17, TIMP3)
concord <- merge(
  cor_df_255[, c("GLUT", "Reg", "r", "p")],
  cor_df_249[, c("GLUT", "Reg", "r", "p")],
  by = c("GLUT", "Reg"), suffixes = c("_255", "_249")
)
concord$same_sign     <- sign(concord$r_255) == sign(concord$r_249)
concord$both_signif   <- concord$p_255 < 0.05 & concord$p_249 < 0.05
concord$agree_signif  <- concord$same_sign & concord$both_signif
write.csv(concord, file.path(outdir_249, "B10_cross_dataset_concordance.csv"),
          row.names = FALSE)
cat("Pairs agreeing in sign:", sum(concord$same_sign), "/",  nrow(concord), "\n")
cat("Pairs significant in both:", sum(concord$both_signif), "\n")




# ─── END ──────────────────────────────────────────────────────────────────────

sink(file.path(outdir_249, "B_sessionInfo.txt"))
sessionInfo(); sink()