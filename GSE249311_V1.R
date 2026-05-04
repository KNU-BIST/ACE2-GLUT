############################################################
# Placental RNA-seq – GSE249311 (cleaned, no Unknown samples)
# Goal:
# 1. Load count matrix
# 2. Create sample metadata
# 3. DESeq2 normalization (VST)
# 4. Extract GLUT + RAS genes
# 5. Heatmap + correlation plot (GLUT vs ACE2/RAS)
############################################################

## 1. Libraries ------------------------------------------------------------
library(readxl)
library(DESeq2)
library(dplyr)
library(pheatmap)
library(ggplot2)

## 2. Load count matrix from Excel ----------------------------------------
counts_df <- read_excel("GSE249311_count_matrix.xlsx")

# first column = gene symbols
gene_ids <- counts_df[[1]]
counts   <- as.matrix(counts_df[ , -1])

rownames(counts) <- make.unique(gene_ids)
mode(counts) <- "numeric"

# sample names from columns
sample_names <- colnames(counts)

## 3. Define sample metadata ----------------------------------------------


condition <- c(rep("Control", 11), rep("GDM_diet", 5), rep("GDM_insulin", 9), rep("T2DM", 5))
# quick check: length must equal number of samples
stopifnot(length(condition) == length(sample_names))

sample_info <- data.frame(
  row.names = sample_names,
  condition = factor(condition)
)

table(sample_info$condition)

## 4. Build DESeq2 object and filter genes --------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = sample_info,
  design    = ~ condition
)

# remove very lowly expressed genes
dds <- dds[ rowSums(counts(dds)) > 10, ]

## 5. Variance-stabilizing transformation ---------------------------------
vsd <- vst(dds)
expr_matrix <- assay(vsd)   # genes x samples

## 6. Define genes of interest (GLUT + RAS panel) -------------------------
glut_genes <- c("SLC2A1","SLC2A3","SLC2A4","SLC2A9")

ras_genes <- c(
  "ACE","ACE2","REN","AGT",
  "AGTR1","AGTR2","MAS1",
  "MME","THOP1","ANPEP","LNPEP",
  "ADAM17","TIMP3"
)

genes_of_interest <- unique(c(glut_genes, ras_genes))

# keep only genes actually present
present_genes <- intersect(genes_of_interest, rownames(expr_matrix))
expr_subset   <- expr_matrix[present_genes, ]

## 7. Order samples by condition for plotting -----------------------------
sample_info$condition <- factor(sample_info$condition)
expr_subset <- expr_subset[ , rownames(sample_info) ]

## 8. Heatmap of GLUT + RAS expression ------------------------------------
annotation_col <- data.frame(Condition = sample_info$condition)
rownames(annotation_col) <- rownames(sample_info)

pheatmap(
  expr_subset,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = annotation_col,
  fontsize_row = 9,
  main = "GLUT and RAS gene expression (GSE249311)"
)

# save high-res version
png("GSE249311_GLUT_RAS_heatmap.png", width = 2000, height = 2000, res = 300)
pheatmap(
  expr_subset,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = annotation_col,
  fontsize_row = 9,
  main = "GLUT and RAS gene expression (GSE249311)"
)
dev.off()
dev.new()

## 9. Correlation: GLUT vs ACE2/RAS ---------------------------------------
present_glut <- intersect(glut_genes, rownames(expr_subset))
present_ras  <- setdiff(rownames(expr_subset), present_glut)

expr_glut <- t(expr_subset[present_glut, ])   # samples x GLUTs
expr_ras  <- t(expr_subset[present_ras, ])    # samples x RAS genes

# Pearson correlation across all samples
cor_mat <- cor(expr_glut, expr_ras,
               use = "pairwise.complete.obs",
               method = "pearson")

# heatmap of correlations
pheatmap(
  cor_mat,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = colorRampPalette(c("blue","white","red"))(100),
  main  = "Correlation: GLUT vs ACE2/RAS genes",
  angle_col = 45
)

png("GSE249311_GLUT_RAS_correlation.png", width = 2000, height = 1800, res = 300)
pheatmap(
  cor_mat,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = colorRampPalette(c("blue","white","red"))(100),
  main  = "Correlation: GLUT vs ACE2/RAS genes",
  angle_col = 45
)
dev.off()

## 10. (Optional) export correlation matrix -------------------------------
write.csv(cor_mat, file = "GSE249311_GLUT_RAS_correlation_matrix.csv")



############################################################
# 11. Structured heatmap: GLUT vs RAS vs RAS regulators
############################################################

# Define groups of genes

# GLUT transporters
glut_genes <- c("SLC2A1","SLC2A3","SLC2A4","SLC2A9")

# Core RAS
ras_core <- c("ACE","ACE2","REN","AGT","AGTR1","AGTR2","MAS1")

# RAS regulators
ras_reg <- c("MME","THOP1","ANPEP","LNPEP","ADAM17","TIMP3")

expr_panel <- expr_matrix[intersect(gene_panel, rownames(expr_matrix)), ]





expr_panel_scaled <- t(scale(t(expr_panel)))

# remove rows with NA
expr_panel_scaled <- expr_panel_scaled[complete.cases(expr_panel_scaled), ]



gene_category <- data.frame(
  Category = c(
    rep("GLUT", length(glut_genes)),
    rep("RAS_core", length(ras_core)),
    rep("RAS_regulator", length(ras_reg))
  )
)

rownames(gene_category) <- gene_panel

# keep only genes present in dataset
gene_category <- gene_category[rownames(expr_panel_scaled), , drop = FALSE]



ann_colors <- list(
  Category = c(
    GLUT = "#1b9e77",
    RAS_core = "#d95f02",
    RAS_regulator = "#7570b3"
  )
)


gaps_row <- c(
  length(glut_genes),
  length(glut_genes) + length(ras_core)
)


pheatmap(
  expr_panel_scaled,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_row = gene_category,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  gaps_row = gaps_row,
  color = colorRampPalette(c("navy","white","firebrick3"))(100),
  fontsize_row = 10,
  main = "Placental RAS–GLUT metabolic axis (GSE249311)"
)






############################################################
# Focused correlations: ACE2, ACE, AGTR1 vs GLUTs
############################################################

library(reshape2)
library(ggplot2)

# 1. Define genes
glut_genes <- c("SLC2A1","SLC2A3","SLC2A4","SLC2A9")
reg_genes  <- c("ACE","ACE2","REN","AGT",
                "AGTR1","AGTR2","MAS1",
                "MME","THOP1","ANPEP","LNPEP",
                "ADAM17","TIMP3")   # you can add/remove here

genes_needed <- c(glut_genes, reg_genes)
present_genes <- intersect(genes_needed, rownames(expr_matrix))
present_genes
# make sure all are present; if not, drop missing ones or rename

# 2. Subset expression matrix (genes x samples)
expr_sub <- expr_matrix[present_genes, ]

# split into GLUT and regulator matrices (samples x genes)
expr_glut <- t(expr_sub[intersect(glut_genes, rownames(expr_sub)), ])
expr_reg  <- t(expr_sub[intersect(reg_genes,  rownames(expr_sub)), ])

# 3. Compute Pearson r and p‑values for each pair

pairs <- expand.grid(
  GLUT = colnames(expr_glut),
  Reg  = colnames(expr_reg),
  stringsAsFactors = FALSE
)

get_cor <- function(g, r) {
  x <- expr_glut[, g]
  y <- expr_reg[, r]
  ct <- cor.test(x, y, method = "pearson")
  data.frame(
    GLUT = g,
    Reg  = r,
    r    = unname(ct$estimate),
    p    = ct$p.value,
    stringsAsFactors = FALSE
  )
}

# use lapply + do.call instead of apply + rbind on a matrix
cor_list <- lapply(seq_len(nrow(pairs)), function(i) {
  get_cor(pairs$GLUT[i], pairs$Reg[i])
})
cor_df <- do.call(rbind, cor_list)

# set factor order
cor_df$GLUT <- factor(cor_df$GLUT, levels = glut_genes)
cor_df$Reg  <- factor(cor_df$Reg,  levels = reg_genes)

# add labels and significance
cor_df$label <- sprintf("%.2f", cor_df$r)
cor_df$sign  <- cut(
  cor_df$p,
  breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
  labels = c("***","**","*","")
)




############################################################
# 4. Heatmap‑style dot plot with r and significance
############################################################

ggplot(cor_df, aes(x = Reg, y = GLUT, fill = r)) +
  geom_tile(color = "grey80") +
  geom_text(aes(label = paste0(label, sign)), size = 4) +
  scale_fill_gradient2(limit = c(-1,1),
                       low = "blue", mid = "white", high = "red",
                       name = "Pearson r") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(
    title = "Correlation between GLUTs and ACE2 / RAS genes",
    x = "RAS genes",
    y = "GLUT genes"
  )

# save high‑resolution version
ggsave("Focused_GLUT_ACE2_correlation.png",
       width = 12, height = 4, dpi = 300)












## 11. GLUT Expression Boxplots -------------------------------------------

shape_values <- c(
  Control = 16,
  GDM_diet = 17,
  GDM_insulin = 15,
  T2DM = 18
)

for(g in glut_present){

  gene_df <- glut_long[glut_long$variable == g, ]

  p <- ggplot(gene_df,
              aes(x = Condition,
                  y = value,
                  shape = Condition)) +

    geom_boxplot(fill = "grey90", color = "black", outlier.shape = NA) +

    geom_jitter(width = 0.15, size = 3) +

    scale_shape_manual(values = shape_values) +

    theme_bw() +

    labs(
      title = paste("Placental", g, "expression"),
      y = "VST normalized expression",
      x = "Condition"
    ) +

    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  ggsave(
    paste0("GLUT_", g, "_expression.png"),
    p,
    width = 5,
    height = 5,
    dpi = 300
  )
}

# save
ggsave(
  "GSE249311_GLUT_expression_boxplots.png",
  glut_plot,
  width = 8,
  height = 6,
  dpi = 300
)


#########Another option###########################
shape_values <- c(
  Control = 16,
  GDM_diet = 17,
  GDM_insulin = 15,
  T2DM = 18
)

# ==== CUSTOM SETTINGS ====
y_text_size <- 14
x_text_size <- 14
title_size  <- 16
point_size  <- 3
box_line    <- 0.8

for(g in glut_present){
  
  gene_df <- glut_long[glut_long$variable == g, ]
  
  p <- ggplot(gene_df,
              aes(x = Condition,
                  y = value,
                  shape = Condition)) +
    
    geom_boxplot(
      fill = "white",
      color = "black",
      linewidth = box_line,
      outlier.shape = NA
    ) +
    
    geom_jitter(
      width = 0.15,
      size = point_size
    ) +
    
    scale_shape_manual(values = shape_values) +
    
    theme_bw() +
    
    labs(
      title = g,   # cleaner than paste()
      y = "VST normalized expression",
      x = NULL     # 🔥 removes x-axis title
    ) +
    
    theme(
      legend.position = "none",
      
      # 🔥 CENTER TITLE
      plot.title = element_text(
        hjust = 0.5,
        size = title_size,
        face = "bold"
      ),
      
      # 🔥 Y-AXIS TEXT BIGGER
      axis.text.y = element_text(size = y_text_size),
      axis.title.y = element_text(size = y_text_size + 1),
      
      # X-axis formatting
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        size = x_text_size
      ),
      
      # Optional: cleaner panel
      panel.grid = element_blank()
    )
  
  ggsave(
    paste0("GLUT_", g, "_expression.png"),
    p,
    width = 5,
    height = 5,
    dpi = 300
  )
}







###############Save the statistics##########################################
library(dplyr)

pairwise_results <- glut_long %>%
  group_by(variable) %>%
  do({
    
    df <- .
    
    comparisons <- list(
      c("Control","GDM_diet"),
      c("Control","GDM_insulin"),
      c("Control","T2DM")
    )
    
    res_list <- lapply(comparisons, function(comp){
      
      sub_df <- df[df$Condition %in% comp, ]
      
      test <- t.test(value ~ Condition, data = sub_df)
      
      data.frame(
        comparison = paste(comp, collapse = "_vs_"),
        p_value = test$p.value
      )
    })
    
    do.call(rbind, res_list)
  })




pairwise_results$FDR <- p.adjust(pairwise_results$p_value, method = "fdr")

pairwise_results

library(tidyr)

final_table <- pairwise_results %>%
  pivot_wider(
    names_from = comparison,
    values_from = c(p_value, FDR)
  )

final_table




########################Final excell###############################
library(openxlsx)

write.xlsx(final_table, "GLUT_pairwise_Control_comparisons.xlsx")






## 12. Differential expression --------------------------------------------

dds <- DESeq(dds)

# example comparison
res_gdm <- results(dds, contrast = c("condition","GDM_insulin","Control"))

# convert to dataframe
res_df <- as.data.frame(res_gdm)

# remove NA
res_df <- na.omit(res_df)

# select significant genes
sig_genes <- rownames(res_df[res_df$padj < 0.05, ])

write.csv(res_df, "GSE249311_DEG_GDMinsulin_vs_Control.csv")

















## 13. Pathway enrichment --------------------------------------------------

library(clusterProfiler)
library(org.Hs.eg.db)

# convert gene symbols to Entrez IDs
gene_entrez <- bitr(
  sig_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# KEGG enrichment
kegg_res <- enrichKEGG(
  gene = gene_entrez$ENTREZID,
  organism = "hsa"
)

# plot enrichment
kegg_plot <- dotplot(kegg_res, showCategory = 15) +
  ggtitle("Pathway enrichment: GDM vs Control placenta")

print(kegg_plot)

# save
ggsave(
  "GSE249311_KEGG_pathway_enrichment.png",
  kegg_plot,
  width = 8,
  height = 6,
  dpi = 300
)