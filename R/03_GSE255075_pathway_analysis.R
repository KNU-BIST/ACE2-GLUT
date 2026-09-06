# =============================================================================
# Phase 3 — Pathway and enrichment analysis of GSE255075
# Script: 03_GSE255075_pathway_analysis.R
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
# Reproduce pathway-level analyses derived from the GSE255075 DESeq2 results,
# including GAGE/KEGG analysis, Pathview pathway rendering, Gene Ontology
# over-representation analysis, GO GSEA, and enrichment-curve generation.
#
# Reproducibility boundary
# ------------------------
# This phase reads persisted differential-expression results from Phase 1.
# It does not refit the DESeq2 model.
#
# Pathway analyses are transcriptomic context and do not establish aspirin,
# insulin, ACE2, or GLUT causal mechanisms in the institutional placental
# cohort.
#
# Inputs
# ------
# outputs/tables/GSE255075/DESeq2_GDM_vs_Normal_shrunk.csv
#
# Outputs
# -------
# outputs/tables/GSE255075/KEGG_*
# outputs/tables/GSE255075/GO_*
# outputs/tables/GSE255075/GSEA_*
# figures/exploratory/GSE255075/GSEA_*
# figures/exploratory/GSE255075/pathview/
# =============================================================================

if (!requireNamespace("here", quietly = TRUE)) {
  stop("Install package 'here' before running this script.")
}
source(here::here("R", "00_project_setup.R"))

de_file <- file.path(
  path_tables_255,
  "DESeq2_GDM_vs_Normal_shrunk.csv"
)

assert_object_exists(
  de_file,
  "GSE255075 differential-expression result"
)

res_df_255 <- read.csv(
  de_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

message("\n===== GSE255075: pathway analysis =====\n")

# -----------------------------------------------------------------------------
# 1. Map fold changes to Entrez IDs
# -----------------------------------------------------------------------------

gene_ids_255 <- res_df_255$gene
is_ensembl_255 <- mean(
  grepl("^ENSG", gene_ids_255)
) > 0.5

fc_df_255 <- res_df_255 |>
  dplyr::select(gene, log2FoldChange) |>
  dplyr::filter(!is.na(log2FoldChange))

if (is_ensembl_255) {
  fc_df_255$gene <- sub(
    "\\..*$",
    "",
    fc_df_255$gene
  )

  map_255 <- clusterProfiler::bitr(
    fc_df_255$gene,
    fromType = "ENSEMBL",
    toType = c("ENTREZID", "SYMBOL"),
    OrgDb = org.Hs.eg.db
  )

  fc_mapped_255 <- fc_df_255 |>
    dplyr::inner_join(
      map_255,
      by = c("gene" = "ENSEMBL")
    ) |>
    dplyr::distinct(ENTREZID, .keep_all = TRUE)
} else {
  map_255 <- clusterProfiler::bitr(
    fc_df_255$gene,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )

  fc_mapped_255 <- fc_df_255 |>
    dplyr::inner_join(
      map_255,
      by = c("gene" = "SYMBOL")
    ) |>
    dplyr::distinct(ENTREZID, .keep_all = TRUE)
}

foldchanges_255 <- setNames(
  fc_mapped_255$log2FoldChange,
  fc_mapped_255$ENTREZID
)

foldchanges_255 <- sort(
  foldchanges_255,
  decreasing = TRUE
)

write.csv(
  fc_mapped_255,
  file.path(
    path_tables_255,
    "KEGG_foldchange_ENTREZ_mapping.csv"
  ),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 2. GAGE KEGG enrichment
# -----------------------------------------------------------------------------

kegg_hsa_255 <- gage::kegg.gsets(species = "hsa")
kegg_sets_hsa <- kegg_hsa_255$kg.sets[
  kegg_hsa_255$sigmet.idx
]

kegg_res_255 <- gage::gage(
  foldchanges_255,
  gsets = kegg_sets_hsa,
  same.dir = TRUE
)

if (!is.null(kegg_res_255$greater)) {
  write.csv(
    as.data.frame(kegg_res_255$greater) |>
      tibble::rownames_to_column("pathway") |>
      dplyr::arrange(p.val),
    file.path(
      path_tables_255,
      "KEGG_GAGE_upregulated.csv"
    ),
    row.names = FALSE
  )
}

if (!is.null(kegg_res_255$less)) {
  write.csv(
    as.data.frame(kegg_res_255$less) |>
      tibble::rownames_to_column("pathway") |>
      dplyr::arrange(p.val),
    file.path(
      path_tables_255,
      "KEGG_GAGE_downregulated.csv"
    ),
    row.names = FALSE
  )
}

# -----------------------------------------------------------------------------
# 3. Pathway membership overlap and Pathview
# -----------------------------------------------------------------------------

goi_map_255 <- clusterProfiler::bitr(
  genes_of_interest,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

goi_entrez_255 <- unique(goi_map_255$ENTREZID)
entrez2sym_255 <- setNames(
  goi_map_255$SYMBOL,
  goi_map_255$ENTREZID
)

hits_255_all <- lapply(
  names(kegg_sets_hsa),
  function(pathway_name) {
    overlap <- intersect(
      goi_entrez_255,
      kegg_sets_hsa[[pathway_name]]
    )

    if (length(overlap) > 0) {
      overlap
    } else {
      NULL
    }
  }
)

names(hits_255_all) <- names(kegg_sets_hsa)
hits_255 <- Filter(
  Negate(is.null),
  hits_255_all
)

hit_tbl_255 <- data.frame(
  pathway = names(hits_255),
  n_genes_present = vapply(
    hits_255,
    length,
    integer(1)
  ),
  genes_present = vapply(
    hits_255,
    function(x) {
      paste(
        entrez2sym_255[x],
        collapse = ","
      )
    },
    character(1)
  )
) |>
  dplyr::arrange(
    dplyr::desc(n_genes_present),
    pathway
  )

write.csv(
  hit_tbl_255,
  file.path(
    path_tables_255,
    "KEGG_pathways_containing_GOI.csv"
  ),
  row.names = FALSE
)

up_paths_255 <- if (!is.null(kegg_res_255$greater)) {
  rownames(kegg_res_255$greater)
} else {
  character(0)
}

down_paths_255 <- if (!is.null(kegg_res_255$less)) {
  rownames(kegg_res_255$less)
} else {
  character(0)
}

path_ids_255 <- unique(
  c(
    substr(head(up_paths_255, 20), 1, 8),
    substr(
      hit_tbl_255$pathway[
        hit_tbl_255$pathway %in% up_paths_255
      ],
      1,
      8
    ),
    substr(
      hit_tbl_255$pathway[
        hit_tbl_255$pathway %in% down_paths_255
      ],
      1,
      8
    )
  )
)

path_ids_255 <- path_ids_255[
  nzchar(path_ids_255)
]

pathview_dir <- file.path(
  path_fig_255,
  "pathview"
)

dir.create(
  pathview_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

old_wd <- getwd()

tryCatch(
  {
    setwd(pathview_dir)

    invisible(
      lapply(
        path_ids_255,
        function(pid) {
          try(
            pathview::pathview(
              gene.data = foldchanges_255,
              pathway.id = pid,
              species = "hsa",
              kegg.native = TRUE,
              same.layer = FALSE,
              out.suffix = "GDM_vs_Normal"
            ),
            silent = TRUE
          )
        }
      )
    )
  },
  finally = {
    setwd(old_wd)
  }
)

writeLines(
  path_ids_255,
  file.path(
    path_tables_255,
    "Pathview_plotted_pathway_ids.txt"
  )
)

# -----------------------------------------------------------------------------
# 4. GO mapping
# -----------------------------------------------------------------------------

if (is_ensembl_255) {
  res_go_filtered <- res_df_255 |>
    dplyr::filter(!is.na(log2FoldChange)) |>
    dplyr::mutate(
      gene = sub("\\..*$", "", gene)
    )

  map_go_255 <- clusterProfiler::bitr(
    res_go_filtered$gene,
    fromType = "ENSEMBL",
    toType = c("ENTREZID", "SYMBOL"),
    OrgDb = org.Hs.eg.db
  )

  res_mapped_go_255 <- res_go_filtered |>
    dplyr::inner_join(
      map_go_255,
      by = c("gene" = "ENSEMBL")
    ) |>
    dplyr::distinct(
      ENTREZID,
      .keep_all = TRUE
    )
} else {
  res_go_filtered <- res_df_255 |>
    dplyr::filter(!is.na(log2FoldChange))

  map_go_255 <- clusterProfiler::bitr(
    res_go_filtered$gene,
    fromType = "SYMBOL",
    toType = c("ENTREZID", "SYMBOL"),
    OrgDb = org.Hs.eg.db
  )

  res_mapped_go_255 <- res_go_filtered |>
    dplyr::inner_join(
      map_go_255,
      by = c("gene" = "SYMBOL")
    ) |>
    dplyr::distinct(
      ENTREZID,
      .keep_all = TRUE
    )
}

write.csv(
  res_mapped_go_255,
  file.path(
    path_tables_255,
    "GO_res_mapped_to_ENTREZ.csv"
  ),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 5. GO over-representation analysis
# -----------------------------------------------------------------------------

padj_go <- 0.05
lfc_go <- 1

deg_all_255 <- res_mapped_go_255 |>
  dplyr::filter(
    !is.na(padj),
    padj < padj_go
  )

deg_up_255 <- deg_all_255 |>
  dplyr::filter(
    log2FoldChange >= lfc_go
  )

deg_down_255 <- deg_all_255 |>
  dplyr::filter(
    log2FoldChange <= -lfc_go
  )

run_enrich_go <- function(entrez_vec, ontology) {
  clusterProfiler::enrichGO(
    gene = entrez_vec,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = ontology,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
  )
}

ego_bp <- run_enrich_go(
  deg_all_255$ENTREZID,
  "BP"
)

ego_mf <- run_enrich_go(
  deg_all_255$ENTREZID,
  "MF"
)

ego_cc <- run_enrich_go(
  deg_all_255$ENTREZID,
  "CC"
)

ego_up <- run_enrich_go(
  deg_up_255$ENTREZID,
  "BP"
)

ego_down <- run_enrich_go(
  deg_down_255$ENTREZID,
  "BP"
)

go_objects <- list(
  BP_all = ego_bp,
  MF_all = ego_mf,
  CC_all = ego_cc,
  BP_up = ego_up,
  BP_down = ego_down
)

for (nm in names(go_objects)) {
  write.csv(
    as.data.frame(go_objects[[nm]]),
    file.path(
      path_tables_255,
      paste0("GO_ORA_", nm, ".csv")
    ),
    row.names = FALSE
  )
}

# -----------------------------------------------------------------------------
# 6. GO GSEA
# -----------------------------------------------------------------------------

gene_list_255 <- sort(
  setNames(
    res_mapped_go_255$log2FoldChange,
    res_mapped_go_255$ENTREZID
  ),
  decreasing = TRUE
)

gsea_bp_255 <- clusterProfiler::gseGO(
  geneList = gene_list_255,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = FALSE
)

gsea_res_255 <- as.data.frame(gsea_bp_255)

write.csv(
  gsea_res_255,
  file.path(
    path_tables_255,
    "GSEA_GO_BP.csv"
  ),
  row.names = FALSE
)

if (nrow(gsea_res_255) > 0) {
  gsea_res_255$Direction <- ifelse(
    gsea_res_255$NES >= 0,
    "Upregulated",
    "Downregulated"
  )

  gsea_bp_plot <- gsea_bp_255
  gsea_bp_plot@result <- gsea_res_255

  p_gsea_split <- enrichplot::dotplot(
    gsea_bp_plot,
    showCategory = 20,
    x = "NES",
    label_format = 45
  ) +
    ggplot2::facet_grid(
      Direction ~ .,
      scales = "free_y"
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::labs(
      title = "GO BP GSEA — GSE255075",
      x = "NES",
      y = NULL
    )

  save_plot(
    p_gsea_split,
    file.path(
      path_fig_255,
      "GSE255075_GSEA_BP_up_vs_down.png"
    ),
    width = 17 / 2.54,
    height = 18 / 2.54
  )

  ggplot2::ggsave(
    file.path(
      path_fig_255,
      "GSE255075_GSEA_BP_up_vs_down.pdf"
    ),
    p_gsea_split,
    width = 17 / 2.54,
    height = 18 / 2.54
  )

  # Per-term enrichment curves for positively enriched terms with padj <= 0.05.
  up_terms_255 <- gsea_res_255 |>
    dplyr::filter(
      !is.na(NES),
      NES > 0,
      !is.na(p.adjust),
      p.adjust <= 0.05
    ) |>
    dplyr::arrange(p.adjust)

  gsea_curve_dir <- file.path(
    path_fig_255,
    "GSEA_curves"
  )

  dir.create(
    gsea_curve_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  for (i in seq_len(nrow(up_terms_255))) {
    term_id <- up_terms_255$ID[i]
    description <- up_terms_255$Description[i]
    nes_value <- round(up_terms_255$NES[i], 2)
    padj_value <- signif(
      up_terms_255$p.adjust[i],
      3
    )

    title_text <- paste0(
      description,
      "\nNES = ",
      nes_value,
      " | adjusted P = ",
      padj_value
    )

    p_enrichment <- tryCatch(
      enrichplot::gseaplot2(
        gsea_bp_255,
        geneSetID = term_id,
        title = title_text,
        pvalue_table = FALSE
      ),
      error = function(e) NULL
    )

    if (is.null(p_enrichment)) {
      next
    }

    safe_id <- gsub(
      "[^A-Za-z0-9_-]",
      "_",
      term_id
    )

    safe_description <- gsub(
      "[^A-Za-z0-9_-]",
      "_",
      substr(description, 1, 60)
    )

    outfile <- file.path(
      gsea_curve_dir,
      paste0(
        sprintf("%03d", i),
        "_",
        safe_description,
        "_",
        safe_id,
        ".png"
      )
    )

    ggplot2::ggsave(
      outfile,
      plot = p_enrichment,
      width = 7.5,
      height = 5.5,
      dpi = dpi_out
    )
  }
}

write_session_info("03_GSE255075_pathway_analysis.R")

message("\nPhase 3 complete.\n")
