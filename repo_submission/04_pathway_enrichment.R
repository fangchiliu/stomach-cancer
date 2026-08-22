# ==============================================================================
# 04_pathway_enrichment.R — GSEA on the TCGA-STAD tumor vs normal contrast
#
# Ranks every gene by direction x significance, then runs GSEA against Hallmark,
# Reactome, GO:BP and KEGG. Ends with a publication-style EMT enrichment figure.
# ==============================================================================

source("00_setup.R")
theme_set_project()

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(enrichplot)
  library(org.Hs.eg.db)
  library(msigdbr)
  library(tidyverse)
})

set.seed(SEED)

# ---- Ranked gene list ---------------------------------------------------------
res_df         <- read.csv(exp_path("DEA.csv"), row.names = 1)
gene_metadata  <- readRDS(exp_path("gene_metadata.rds"))
res_df$gene_id <- rownames(res_df)
res_df         <- left_join(res_df, gene_metadata, by = "gene_id")

# sign(log2FC) x -log10(p) keeps both direction and strength of evidence, which
# a plain fold-change ranking throws away.
gsea_input <- res_df %>%
  filter(!is.na(padj), !is.na(log2FoldChange), !is.na(gene_name)) %>%
  mutate(rank_metric = sign(log2FoldChange) * -log10(pvalue + 1e-300)) %>%
  arrange(desc(rank_metric))

symbol_to_entrez <- bitr(gsea_input$gene_name,
                         fromType = "SYMBOL", toType = "ENTREZID",
                         OrgDb = org.Hs.eg.db)

gsea_input <- gsea_input %>%
  left_join(symbol_to_entrez, by = c("gene_name" = "SYMBOL")) %>%
  filter(!is.na(ENTREZID)) %>%
  distinct(ENTREZID, .keep_all = TRUE)

gene_list <- setNames(gsea_input$rank_metric, gsea_input$ENTREZID)
saveRDS(gsea_input, exp_path("gsea_input.rds"))
cat("Genes in ranked list:", length(gene_list), "\n")

# ---- MSigDB collections -------------------------------------------------------
get_sets <- function(category, subcategory = NULL) {
  msigdbr(species = "Homo sapiens",
          category = category, subcategory = subcategory) %>%
    dplyr::select(gs_name, entrez_gene) %>%
    mutate(entrez_gene = as.character(entrez_gene))
}

h_sets  <- get_sets("H")
c2_sets <- get_sets("C2", "CP:REACTOME")
c5_sets <- get_sets("C5", "GO:BP")

# ---- Run GSEA ------------------------------------------------------------------
run_gsea <- function(ranks, sets, label) {
  message("Running GSEA: ", label)
  GSEA(geneList = ranks, TERM2GENE = sets,
       pvalueCutoff = PADJ_CUTOFF, pAdjustMethod = "BH",
       minGSSize = 15, maxGSSize = 500,
       eps = 0, seed = TRUE, verbose = FALSE)
}

gsea_hallmark <- cache_rds("gsea_hallmark.rds", run_gsea(gene_list, h_sets,  "Hallmark"))
gsea_reactome <- cache_rds("gsea_reactome.rds", run_gsea(gene_list, c2_sets, "Reactome"))
gsea_gobp     <- cache_rds("gsea_gobp.rds",     run_gsea(gene_list, c5_sets, "GO:BP"))

gsea_kegg <- cache_rds("gsea_kegg.rds", gseKEGG(
  geneList = gene_list, organism = "hsa",
  minGSSize = 15, maxGSSize = 500,
  pvalueCutoff = PADJ_CUTOFF, pAdjustMethod = "BH",
  eps = 0, verbose = FALSE
))

for (nm in c("hallmark", "reactome", "gobp", "kegg")) {
  obj <- get(paste0("gsea_", nm))
  write.csv(as.data.frame(obj),
            exp_path(sprintf("GSEA_%s.csv", toupper(nm))), row.names = FALSE)
  cat(sprintf("%-9s significant pathways: %d\n", nm, nrow(as.data.frame(obj))))
}

# ---- Dotplots ------------------------------------------------------------------
save_dotplot <- function(gsea_obj, title, filename, n = 20, text_size = 8) {
  png(plot_path(filename), width = 2800, height = 2200, res = 300)
  print(
    dotplot(gsea_obj, showCategory = n, split = ".sign") +
      facet_grid(. ~ .sign) +
      ggtitle(title) +
      theme(axis.text.y = element_text(size = text_size))
  )
  dev.off()
}

save_dotplot(gsea_hallmark, "GSEA — Hallmark (Tumor vs Normal)",  "GSEA_hallmark_dotplot.png")
save_dotplot(gsea_reactome, "GSEA — Reactome (Tumor vs Normal)",  "GSEA_reactome_dotplot.png", text_size = 7)
save_dotplot(gsea_kegg,     "GSEA — KEGG (Tumor vs Normal)",      "GSEA_KEGG_dotplot.png")

# ---- NES bar charts -------------------------------------------------------------
nes_barplot <- function(gsea_obj, label_col, strip_prefix = NULL, title) {
  df <- as.data.frame(gsea_obj) %>%
    arrange(NES) %>%
    mutate(
      label = .data[[label_col]],
      label = if (!is.null(strip_prefix)) str_remove(label, strip_prefix) else label,
      label = str_replace_all(label, "_", " "),
      Direction = ifelse(NES > 0, "Activated in Tumor", "Suppressed in Tumor")
    )

  ggplot(df, aes(x = reorder(label, NES), y = NES, fill = Direction)) +
    geom_col() +
    scale_fill_manual(values = c("Activated in Tumor"  = COL_TUMOR,
                                 "Suppressed in Tumor" = COL_NORMAL)) +
    coord_flip() +
    labs(title = title, subtitle = "Gastric cancer tumor vs normal tissue",
         x = NULL, y = "Normalized Enrichment Score (NES)", fill = NULL) +
    theme(legend.position = "bottom")
}

save_plot(nes_barplot(gsea_hallmark, "ID", "HALLMARK_",
                      "GSEA Hallmark — Normalized Enrichment Scores"),
          "GSEA_hallmark_NES_barplot.png", width = 10, height = 8)

save_plot(nes_barplot(gsea_kegg, "Description", NULL,
                      "GSEA KEGG — Normalized Enrichment Scores"),
          "GSEA_KEGG_NES_barplot.png", width = 10, height = 8)

# ---- Hallmark network and ridge views --------------------------------------------
png(plot_path("GSEA_hallmark_emap.png"), width = 2800, height = 2400, res = 300)
print(emapplot(pairwise_termsim(gsea_hallmark), showCategory = 25) +
        ggtitle("Hallmark pathway similarity network"))
dev.off()

png(plot_path("GSEA_hallmark_ridge.png"), width = 2600, height = 2400, res = 300)
print(ridgeplot(gsea_hallmark, showCategory = 20, fill = "p.adjust") +
        ggtitle("Hallmark — ranked gene distribution per pathway") +
        theme(axis.text.y = element_text(size = 8)))
dev.off()

# ---- Enrichment curves for the strongest hits in each direction -------------------
top_ids <- function(gsea_obj, n = 3) {
  df <- as.data.frame(gsea_obj)
  c(df %>% filter(NES > 0) %>% arrange(p.adjust) %>% slice(1:n) %>% pull(ID),
    df %>% filter(NES < 0) %>% arrange(p.adjust) %>% slice(1:n) %>% pull(ID))
}

png(plot_path("GSEA_hallmark_enrichment_curves.png"),
    width = 3600, height = 2400, res = 300)
print(gseaplot2(gsea_hallmark, geneSetID = top_ids(gsea_hallmark),
                title = "Top activated & suppressed Hallmark pathways",
                pvalue_table = TRUE))
dev.off()

# ---- KEGG pathway diagrams with fold changes overlaid ------------------------------
if (requireNamespace("pathview", quietly = TRUE)) {
  library(pathview)
  fc_vector <- setNames(gsea_input$log2FoldChange, gsea_input$ENTREZID)

  for (pid in top_ids(gsea_kegg)) {
    tryCatch(
      pathview(gene.data = fc_vector, pathway.id = pid, species = "hsa",
               out.suffix = "STAD", kegg.dir = plots_dir),
      error = function(e) message("pathview failed for ", pid, ": ", e$message)
    )
  }
}

# ==============================================================================
# EMT figure — the headline pathway result
# ==============================================================================
emt_id <- "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
stopifnot(emt_id %in% gsea_hallmark@result$ID)

emt_stats <- gsea_hallmark@result[gsea_hallmark@result$ID == emt_id, ]
print(emt_stats[, c("ID", "NES", "pvalue", "p.adjust", "setSize")])

annotation <- paste(
  sprintf("NES = %.2f",  emt_stats$NES),
  sprintf("FDR = %.2e",  emt_stats$p.adjust),
  sprintf("%d genes",    emt_stats$setSize),
  sep = "\n"
)

# Shorten the description that gseaplot2 bakes into the panel.
gsea_plot_obj <- gsea_hallmark
gsea_plot_obj@result[emt_id, "Description"] <- "Epithelial-Mesenchymal Transition"

emt_plot <- gseaplot2(
  gsea_plot_obj, geneSetID = emt_id, title = NULL, color = COL_TUMOR,
  base_size = 14, rel_heights = c(1.6, 0.35, 0.9),
  subplots = 1:3, pvalue_table = FALSE, ES_geom = "line"
)

# Panel 1: running enrichment score
emt_plot[[1]] <- emt_plot[[1]] +
  annotate("text", x = Inf, y = Inf, label = annotation,
           hjust = 1.15, vjust = 1.4, size = 4.6, lineheight = 1.35,
           colour = "#1D3557", fontface = "bold") +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey55", linewidth = 0.4) +
  labs(title    = "Epithelial-Mesenchymal Transition is activated in gastric tumors",
       subtitle = "GSEA · Hallmark gene set · TCGA-STAD tumor vs adjacent normal",
       y        = "Running enrichment score") +
  theme_classic(base_size = 14) +
  theme(
    plot.title    = element_text(face = "bold", size = 16, colour = "#1D3557"),
    plot.subtitle = element_text(size = 11.5, colour = "grey35",
                                 margin = margin(b = 12)),
    axis.line.x   = element_blank(),
    axis.ticks.x  = element_blank(),
    axis.text.x   = element_blank(),
    panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
    legend.position = "none"
  )

# Panel 2: gene hit ticks
emt_plot[[2]] <- emt_plot[[2]] + theme_void()

# Panel 3: ranked metric
emt_plot[[3]] <- emt_plot[[3]] +
  labs(x = "Gene rank (tumor-enriched -> normal-enriched)", y = "Ranked metric") +
  theme_classic(base_size = 14) +
  theme(panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3))

save_plot(emt_plot, "GSEA_EMT_enrichment_curve.png", width = 9, height = 7)
ggsave(plot_path("GSEA_EMT_enrichment_curve.pdf"), emt_plot,
       width = 9, height = 7, bg = "white")

cat("Pathway enrichment complete.\n")
