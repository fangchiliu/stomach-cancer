# ==============================================================================
# 03_differential_expression.R — DESeq2: tumor vs adjacent normal (TCGA-STAD)
#
# Fits the primary differential expression model, exports the result table, and
# saves a VST-normalised matrix that the survival and INHBA scripts reuse.
# ==============================================================================

source("00_setup.R")
theme_set_project()

suppressPackageStartupMessages({
  library(DESeq2)
  library(SummarizedExperiment)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})

rna        <- readRDS(exp_path("stad_rna.rds"))
rna_matrix <- assay(rna, "unstranded")

# ---- Build the DESeq2 object -------------------------------------------------
dds_input <- DESeqDataSetFromMatrix(
  countData = rna_matrix,
  colData   = colData(rna),
  design    = ~ sample_type
)

# Normal as reference, so log2FC reads as tumor over normal.
dds_input$sample_type <- relevel(dds_input$sample_type,
                                 ref = "Solid Tissue Normal")

# Drop genes with almost no support anywhere in the cohort.
dds_input <- dds_input[rowSums(counts(dds_input)) >= 10, ]

cat("Genes tested:", nrow(dds_input), "| samples:", ncol(dds_input), "\n")

# ---- Fit (slow: several minutes on the full cohort) ---------------------------
dds <- cache_rds("stad_dds.rds", DESeq(dds_input))

res    <- results(dds)
res_df <- as.data.frame(res)
summary(res)

write.csv(res_df, exp_path("DEA.csv"))

# ---- VST matrix for downstream per-sample work --------------------------------
# Variance-stabilised values are what the survival split and cross-cohort
# expression plots need; raw counts are not comparable across samples.
vsd <- vst(dds, blind = FALSE)
saveRDS(assay(vsd),                  exp_path("TCGA_vst_matrix.rds"))
saveRDS(as.data.frame(colData(dds)), exp_path("TCGA_sample_metadata.rds"))

# ---- Attach gene symbols ------------------------------------------------------
gene_metadata  <- readRDS(exp_path("gene_metadata.rds"))
res_df$gene_id <- rownames(res_df)
res_df         <- left_join(res_df, gene_metadata, by = "gene_id")

# ---- Volcano plot -------------------------------------------------------------
plot_df <- res_df %>%
  filter(!is.na(padj), !is.na(log2FoldChange), !is.na(gene_name)) %>%
  mutate(
    # Cap the y-axis so a handful of near-zero p-values don't flatten everything.
    neglog10 = pmin(-log10(padj), 100),
    status = case_when(
      padj < PADJ_CUTOFF & log2FoldChange >  LFC_CUTOFF ~ "Up in tumor",
      padj < PADJ_CUTOFF & log2FoldChange < -LFC_CUTOFF ~ "Down in tumor",
      TRUE                                              ~ "Not significant"
    ),
    status = factor(status, levels = c("Up in tumor", "Down in tumor",
                                       "Not significant"))
  )

inhba <- plot_df %>% filter(gene_name == "INHBA")

pal <- c("Up in tumor"     = COL_TUMOR,
         "Down in tumor"   = COL_NORMAL,
         "Not significant" = COL_NEUTRAL)

p_volcano <- ggplot(plot_df, aes(log2FoldChange, neglog10)) +
  geom_point(aes(color = status), size = 1, alpha = 0.55, stroke = 0) +
  geom_vline(xintercept = c(-LFC_CUTOFF, LFC_CUTOFF), linetype = "dashed",
             color = "grey55", linewidth = 0.3) +
  geom_hline(yintercept = -log10(PADJ_CUTOFF), linetype = "dashed",
             color = "grey55", linewidth = 0.3) +
  geom_point(data = inhba, size = 5, shape = 21,
             fill = COL_ACCENT, color = "white", stroke = 1.1) +
  geom_text_repel(
    data = inhba, aes(label = gene_name),
    size = 5, fontface = "bold", color = "#1D1D1D",
    box.padding = 1.2, point.padding = 0.5,
    nudge_x = 1.5, nudge_y = 6,
    segment.color = "#1D1D1D", segment.size = 0.5, min.segment.length = 0
  ) +
  scale_color_manual(values = pal, name = NULL) +
  scale_x_continuous(breaks = seq(-6, 10, 2)) +
  labs(
    title    = "INHBA is significantly upregulated in gastric cancer",
    subtitle = "TCGA-STAD · tumor vs normal · DESeq2",
    x = expression(Log[2]~fold~change),
    y = expression(-Log[10]~adjusted~italic(P))
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  theme(legend.position = "top")

save_plot(p_volcano, "TCGA_volcano_INHBA.png", width = 9, height = 7)

# ---- Significant hits ---------------------------------------------------------
sig <- res_df %>%
  filter(!is.na(padj), padj < PADJ_CUTOFF, abs(log2FoldChange) > LFC_CUTOFF) %>%
  arrange(padj)

write.csv(sig, exp_path("DEA_significant.csv"), row.names = FALSE)
cat("Significant genes (padj <", PADJ_CUTOFF, ", |log2FC| >", LFC_CUTOFF, "):",
    nrow(sig), "\n")
