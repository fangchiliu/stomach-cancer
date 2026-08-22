# ==============================================================================
# 07_inhba_characterization.R — Focused characterization of INHBA
#
# INHBA came out of the two-cohort intersection as the strongest candidate.
# This script (1) confirms upregulation in both cohorts and (2) asks what
# INHBA's expression program is embedded in, via co-expression-driven GSEA.
# ==============================================================================

source("00_setup.R")
theme_set_project()

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
  library(clusterProfiler)
  library(msigdbr)
})

GENE <- "INHBA"
set.seed(SEED)

gene_metadata <- readRDS(exp_path("gene_metadata.rds"))
tcga_vst_mat  <- readRDS(exp_path("TCGA_vst_matrix.rds"))
tcga_coldata  <- readRDS(exp_path("TCGA_sample_metadata.rds"))
expr_df       <- readRDS(exp_path("GSE66229_expr_df.rds"))
pheno         <- readRDS(exp_path("GSE66229_pheno.rds"))

gene_id <- gene_metadata$gene_id[gene_metadata$gene_name == GENE]

# ==============================================================================
# 1. Tumor vs normal in both cohorts
# ==============================================================================
tcga_long <- data.frame(
  expression = as.numeric(tcga_vst_mat[gene_id, ]),
  sample     = colnames(tcga_vst_mat)
) %>%
  mutate(
    group  = ifelse(tcga_coldata[sample, "sample_type"] == "Solid Tissue Normal",
                    "Normal", "Tumor"),
    cohort = "TCGA-STAD (RNA-seq)"
  ) %>%
  dplyr::select(group, cohort, expression)

expr_mat <- expr_df %>% column_to_rownames("gene_symbol") %>% as.matrix()

geo_long <- data.frame(
  expression    = as.numeric(expr_mat[GENE, ]),
  geo_accession = colnames(expr_mat)
) %>%
  left_join(pheno %>% dplyr::select(geo_accession, sample_type), by = "geo_accession") %>%
  mutate(group = sample_type, cohort = "GSE66229 (microarray)") %>%
  dplyr::select(group, cohort, expression)

gene_long <- bind_rows(tcga_long, geo_long) %>%
  filter(!is.na(group)) %>%
  mutate(
    group  = factor(group, levels = c("Normal", "Tumor")),
    cohort = factor(cohort, levels = c("TCGA-STAD (RNA-seq)",
                                       "GSE66229 (microarray)"))
  )

p_expr <- ggplot(gene_long, aes(group, expression, fill = group)) +
  geom_violin(alpha = 0.25, color = NA, width = 0.9, trim = FALSE) +
  geom_boxplot(width = 0.22, outlier.shape = NA, alpha = 0.9,
               color = "grey25", linewidth = 0.4) +
  geom_jitter(aes(color = group), width = 0.12, size = 0.7, alpha = 0.35) +
  stat_compare_means(method = "wilcox.test", label = "p.signif",
                     comparisons = list(c("Normal", "Tumor")),
                     label.y.npc = 0.92, size = 5, bracket.size = 0.4) +
  facet_wrap(~ cohort, scales = "free_y") +
  scale_fill_manual(values  = c("Normal" = COL_NORMAL, "Tumor" = COL_TUMOR)) +
  scale_color_manual(values = c("Normal" = COL_NORMAL, "Tumor" = COL_TUMOR)) +
  labs(
    title    = sprintf("%s is upregulated in gastric cancer across two independent cohorts", GENE),
    subtitle = "Tumor vs adjacent normal · TCGA-STAD (RNA-seq) and GSE66229 (ACRG microarray)",
    caption  = "Wilcoxon rank-sum test:  ns P > 0.05   * P <= 0.05   ** P <= 0.01   *** P <= 0.001   **** P <= 0.0001",
    x = NULL, y = paste(GENE, "expression")
  ) +
  theme(
    legend.position    = "none",
    plot.caption       = element_text(hjust = 0, color = "grey40", size = 10),
    strip.text         = element_text(face = "bold", size = 12),
    panel.grid.major.x = element_blank(),
    panel.spacing      = unit(1.4, "lines")
  )

save_plot(p_expr, sprintf("%s_tumor_vs_normal_both_cohorts.png", GENE),
          width = 10, height = 5.8)

# Effect sizes for the text.
gene_long %>%
  group_by(cohort, group) %>%
  summarise(n = n(), median = median(expression), .groups = "drop") %>%
  print()

# ==============================================================================
# 2. Co-expression GSEA — guilt by association
#
# Rank every gene by its correlation with INHBA across tumors only, then run
# GSEA on that ranking. This asks which programs INHBA moves with inside a
# tumor, which is a different question from what separates tumor from normal.
# ==============================================================================
tumor_samples <- rownames(tcga_coldata)[tcga_coldata$sample_type == "Primary Tumor"]
mat <- tcga_vst_mat[, colnames(tcga_vst_mat) %in% tumor_samples]

# Map Ensembl IDs to symbols, keeping the highest-variance probe per symbol.
sym  <- gene_metadata$gene_name[match(rownames(mat), gene_metadata$gene_id)]
keep <- !is.na(sym) & sym != ""
mat  <- mat[keep, ]
sym  <- sym[keep]

vars <- apply(mat, 1, var)
ord  <- order(sym, -vars)
mat  <- mat[ord, ]
sym  <- sym[ord]
mat  <- mat[!duplicated(sym), ]
rownames(mat) <- sym[!duplicated(sym)]

cat("Tumor expression matrix:", nrow(mat), "genes x", ncol(mat), "samples\n")

# Spearman keeps a handful of extreme samples from driving the ranking.
gene_vec <- mat[GENE, ]
cors <- apply(mat, 1, function(g) cor(g, gene_vec, method = "spearman"))
cors <- sort(cors[names(cors) != GENE], decreasing = TRUE)

cat("Top", GENE, "-correlated genes:\n")
print(head(cors, 15))

h_sets <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)
c2_reactome <- msigdbr(species = "Homo sapiens", category = "C2",
                       subcategory = "CP:REACTOME") %>%
  dplyr::select(gs_name, gene_symbol)

run_gsea <- function(ranks, sets, label) {
  message("Running co-expression GSEA: ", label)
  GSEA(geneList = ranks, TERM2GENE = sets,
       pvalueCutoff = PADJ_CUTOFF, minGSSize = 15, maxGSSize = 500,
       eps = 0, seed = TRUE, verbose = FALSE)
}

gsea_coexpr_h  <- run_gsea(cors, h_sets,      "Hallmark")
gsea_coexpr_re <- run_gsea(cors, c2_reactome, "Reactome")

write.csv(as.data.frame(gsea_coexpr_h),
          exp_path(sprintf("%s_coexpr_GSEA_Hallmark.csv", GENE)), row.names = FALSE)
write.csv(as.data.frame(gsea_coexpr_re),
          exp_path(sprintf("%s_coexpr_GSEA_Reactome.csv", GENE)), row.names = FALSE)

# ---- Pathways positively associated with INHBA ------------------------------------
hd <- as.data.frame(gsea_coexpr_h) %>%
  filter(NES > 0) %>%
  arrange(desc(NES)) %>%
  slice(1:15) %>%
  mutate(pathway = str_to_title(str_replace_all(str_remove(ID, "HALLMARK_"), "_", " ")))

p_pathways <- ggplot(hd, aes(NES, reorder(pathway, NES), size = setSize, color = NES)) +
  geom_point() +
  scale_color_gradient(low = COL_NORMAL, high = COL_TUMOR, name = "NES") +
  scale_size(range = c(3, 9), name = "Gene set size") +
  labs(
    title    = sprintf("Pathways %s is co-regulated with in gastric cancer", GENE),
    subtitle = sprintf("TCGA-STAD tumors · genes ranked by Spearman correlation with %s · Hallmark GSEA", GENE),
    x = "Normalized enrichment score (NES)", y = NULL
  ) +
  theme(panel.grid.major.y = element_blank())

save_plot(p_pathways, sprintf("%s_coexpression_pathways.png", GENE),
          width = 10, height = 7)

cat(GENE, "characterization complete.\n")
