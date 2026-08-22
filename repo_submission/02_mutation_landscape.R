# ==============================================================================
# 02_mutation_landscape.R — Somatic mutation landscape of TCGA-STAD
#
# Summarises the MAF cached by 01_tcga_download.R: cohort-wide mutation burden,
# the most recurrently mutated genes, and the transition/transversion spectrum.
# ==============================================================================

source("00_setup.R")

suppressPackageStartupMessages({
  library(maftools)
})

maf <- readRDS(exp_path("stad_maf.rds"))

# ---- Cohort summary ----------------------------------------------------------
plotmafSummary(maf = maf, rmOutlier = TRUE, dashboard = TRUE)

# ---- Oncoplot: 20 most frequently mutated genes -------------------------------
png(plot_path("oncoplot_top20_mutated_genes.png"),
    width = 2600, height = 2000, res = 300)
oncoplot(maf = maf, top = 20)
dev.off()

# ---- Mutation spectrum -------------------------------------------------------
# Ti/Tv ratio is a standard sanity check on variant calling and hints at the
# dominant mutational process in the cohort.
png(plot_path("titv_spectrum.png"), width = 2600, height = 2000, res = 300)
titv(maf = maf, plot = TRUE, useSyn = TRUE)
dev.off()

# ---- Per-sample burden table --------------------------------------------------
sample_summary <- getSampleSummary(maf)
write.csv(sample_summary, exp_path("STAD_mutation_burden_per_sample.csv"),
          row.names = FALSE)

gene_summary <- getGeneSummary(maf)
write.csv(gene_summary, exp_path("STAD_mutated_genes.csv"), row.names = FALSE)

cat("Median mutations per sample:", median(sample_summary$total), "\n")
cat("Saved mutation summary tables to exports/\n")
