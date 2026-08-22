# ==============================================================================
# 06_validation_gse66229.R — Independent validation in the ACRG cohort
#
# Downloads GSE66229 (ACRG gastric cancer microarray), runs limma tumor vs
# normal, then intersects the hit list with the TCGA DESeq2 results to find
# genes that replicate in both direction and significance across platforms.
# ==============================================================================

source("00_setup.R")
theme_set_project()

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(limma)
  library(tidyverse)
  library(EnhancedVolcano)
})

GEO_ID <- "GSE66229"

# ==============================================================================
# 1. Download and tidy the GEO series
# ==============================================================================
eset <- cache_rds(paste0(GEO_ID, "_raw_eset.rds"), {
  getGEO(GEO_ID, GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]
})

cat("Dimensions:", nrow(eset), "features x", ncol(eset), "samples\n")

pheno <- pData(eset) %>%
  dplyr::select(geo_accession, title, `patient:ch1`, `tissue:ch1`) %>%
  dplyr::rename(patient_id = `patient:ch1`, tissue_type = `tissue:ch1`) %>%
  mutate(
    sample_type = case_when(
      str_detect(tissue_type, regex("tumor",  ignore_case = TRUE)) ~ "Tumor",
      str_detect(tissue_type, regex("normal", ignore_case = TRUE)) ~ "Normal",
      TRUE ~ NA_character_
    )
  )

print(table(pheno$sample_type, useNA = "ifany"))

# ---- Probe to gene mapping ------------------------------------------------------
# Multi-mapping probes list several symbols separated by "///"; keep the first.
probe_map <- fData(eset) %>%
  dplyr::select(ID, `Gene symbol`) %>%
  dplyr::rename(probe_id = ID, gene_symbol = `Gene symbol`) %>%
  mutate(gene_symbol = str_split(gene_symbol, "///", simplify = TRUE)[, 1]) %>%
  filter(gene_symbol != "")

cat("Probes with a symbol:", nrow(probe_map),
    "-> unique genes:", length(unique(probe_map$gene_symbol)), "\n")

# Collapse probes to gene level by averaging.
expr_df <- as.data.frame(exprs(eset)) %>%
  rownames_to_column("probe_id") %>%
  inner_join(probe_map, by = "probe_id") %>%
  dplyr::select(-probe_id) %>%
  group_by(gene_symbol) %>%
  summarise(across(everything(), mean), .groups = "drop")

saveRDS(pheno,   exp_path(paste0(GEO_ID, "_pheno.rds")))
saveRDS(expr_df, exp_path(paste0(GEO_ID, "_expr_df.rds")))

expr_mat <- expr_df %>% column_to_rownames("gene_symbol") %>% as.matrix()
cat("Gene-level matrix:", nrow(expr_mat), "genes x", ncol(expr_mat), "samples\n")

# ==============================================================================
# 2. Differential expression with limma
# ==============================================================================
pheno_aligned <- pheno %>%
  filter(geo_accession %in% colnames(expr_mat)) %>%
  arrange(match(geo_accession, colnames(expr_mat)))

stopifnot(all(pheno_aligned$geo_accession == colnames(expr_mat)))

group  <- factor(pheno_aligned$sample_type, levels = c("Normal", "Tumor"))
design <- model.matrix(~ group)
colnames(design) <- c("Intercept", "Tumor_vs_Normal")

fit <- eBayes(lmFit(expr_mat, design))

res_geo <- topTable(fit, coef = "Tumor_vs_Normal", number = Inf, sort.by = "P") %>%
  rownames_to_column("gene_symbol")

cat("Significant (adj.P < 0.05):", sum(res_geo$adj.P.Val < PADJ_CUTOFF), "\n")
cat("  up in tumor  (logFC >  1):",
    sum(res_geo$adj.P.Val < PADJ_CUTOFF & res_geo$logFC >  1), "\n")
cat("  down in tumor (logFC < -1):",
    sum(res_geo$adj.P.Val < PADJ_CUTOFF & res_geo$logFC < -1), "\n")

write.csv(res_geo, exp_path(paste0(GEO_ID, "_DEA_limma.csv")), row.names = FALSE)

png(plot_path(paste0(GEO_ID, "_volcano.png")), width = 2600, height = 2000, res = 300)
EnhancedVolcano(
  res_geo, lab = res_geo$gene_symbol,
  x = "logFC", y = "adj.P.Val",
  pCutoff = PADJ_CUTOFF, FCcutoff = LFC_ARRAY,
  title = paste0(GEO_ID, ": Tumor vs Normal"),
  subtitle = "ACRG gastric cancer cohort (validation)"
)
dev.off()

# ==============================================================================
# 3. Cross-cohort intersection
# ==============================================================================
res_df         <- read.csv(exp_path("DEA.csv"), row.names = 1)
res_df$gene_id <- rownames(res_df)
gene_metadata  <- readRDS(exp_path("gene_metadata.rds"))
res_df         <- left_join(res_df, gene_metadata, by = "gene_id")

# LFC_ARRAY on both sides: microarray intensities compress fold changes, so a
# common threshold keeps the comparison symmetric.
tcga_sig <- res_df %>%
  filter(!is.na(padj), !is.na(log2FoldChange),
         padj < PADJ_CUTOFF, abs(log2FoldChange) > LFC_ARRAY) %>%
  dplyr::select(gene_name, log2FoldChange, padj) %>%
  dplyr::rename(logFC_TCGA = log2FoldChange, padj_TCGA = padj)

geo_sig <- res_geo %>%
  filter(adj.P.Val < PADJ_CUTOFF, abs(logFC) > LFC_ARRAY) %>%
  dplyr::select(gene_symbol, logFC, adj.P.Val) %>%
  dplyr::rename(gene_name = gene_symbol, logFC_GEO = logFC, padj_GEO = adj.P.Val)

cat("TCGA significant:", nrow(tcga_sig), "| GEO significant:", nrow(geo_sig), "\n")

validated <- inner_join(tcga_sig, geo_sig, by = "gene_name") %>%
  mutate(
    direction_TCGA = ifelse(logFC_TCGA > 0, "Up", "Down"),
    direction_GEO  = ifelse(logFC_GEO  > 0, "Up", "Down"),
    concordant     = direction_TCGA == direction_GEO
  )

validated_concordant <- validated %>%
  filter(concordant) %>%
  mutate(combined_rank = padj_TCGA * padj_GEO) %>%
  arrange(combined_rank)

cat("Significant in both cohorts:", nrow(validated), "\n")
cat("  concordant direction:", sum(validated$concordant), "\n")
cat("  discordant direction:", sum(!validated$concordant), "\n")

write.csv(validated,            exp_path("TCGA_GSE66229_validated_targets.csv"),    row.names = FALSE)
write.csv(validated_concordant, exp_path("TCGA_GSE66229_validated_concordant.csv"), row.names = FALSE)

# ---- Fold-change concordance ------------------------------------------------------
p_corr <- ggplot(validated, aes(logFC_TCGA, logFC_GEO, color = concordant)) +
  geom_point(alpha = 0.7, size = 2.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.6) +
  scale_color_manual(values = c("TRUE" = "#2A9D8F", "FALSE" = "#E76F51"),
                     labels = c("TRUE" = "Concordant", "FALSE" = "Discordant")) +
  labs(
    title    = "Fold change concordance: TCGA-STAD vs GSE66229",
    subtitle = paste0("n = ", nrow(validated), " genes significant in both cohorts"),
    x = "log2FC — TCGA-STAD (RNA-seq)",
    y = "log2FC — GSE66229 (microarray, ACRG cohort)",
    color = NULL
  )

save_plot(p_corr, "TCGA_GSE66229_concordance_scatter.png", width = 8, height = 7)

cat("Pearson r:", round(cor(validated$logFC_TCGA, validated$logFC_GEO), 3), "\n")

# ==============================================================================
# 4. Per-sample expression of the top validated targets
# ==============================================================================
TOP_GENES <- c("GKN1", "GKN2", "ATP4A", "MMP3", "CXCL8", "INHBA")

tcga_vst_mat <- readRDS(exp_path("TCGA_vst_matrix.rds"))
tcga_coldata <- readRDS(exp_path("TCGA_sample_metadata.rds"))

top_gene_ids <- gene_metadata %>%
  filter(gene_name %in% TOP_GENES) %>%
  dplyr::select(gene_id, gene_name)

tcga_expr <- tcga_vst_mat[top_gene_ids$gene_id, , drop = FALSE]
rownames(tcga_expr) <- top_gene_ids$gene_name[
  match(rownames(tcga_expr), top_gene_ids$gene_id)
]

tcga_long <- as.data.frame(tcga_expr) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "sample", values_to = "expression") %>%
  mutate(
    group  = ifelse(tcga_coldata[sample, "sample_type"] == "Solid Tissue Normal",
                    "Normal", "Tumor"),
    cohort = "TCGA-STAD (RNA-seq)"
  ) %>%
  dplyr::select(Gene, group, cohort, expression)

geo_long <- as.data.frame(expr_mat[rownames(expr_mat) %in% TOP_GENES, , drop = FALSE]) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "geo_accession", values_to = "expression") %>%
  left_join(pheno %>% dplyr::select(geo_accession, sample_type), by = "geo_accession") %>%
  mutate(group = sample_type, cohort = "GSE66229 (microarray)") %>%
  dplyr::select(Gene, group, cohort, expression)

combined_long <- bind_rows(tcga_long, geo_long) %>%
  filter(!is.na(group)) %>%
  mutate(Gene  = factor(Gene, levels = TOP_GENES),
         group = factor(group, levels = c("Normal", "Tumor")))

p_box <- ggplot(combined_long, aes(group, expression, fill = group)) +
  geom_boxplot(outlier.size = 0.6, width = 0.6) +
  facet_grid(cohort ~ Gene, scales = "free_y") +
  scale_fill_manual(values = c("Normal" = COL_NORMAL, "Tumor" = COL_TUMOR)) +
  labs(
    title    = "Validated gastric cancer targets: tumor vs normal expression",
    subtitle = "Concordant across TCGA-STAD (RNA-seq) and GSE66229 (ACRG microarray)",
    x = NULL, y = "Expression (VST / log2 intensity)", fill = NULL
  ) +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold", size = 9))

save_plot(p_box, "validated_targets_boxplot.png", width = 14, height = 6)

cat("Validation analysis complete.\n")
