# ==============================================================================
# 01_tcga_download.R — Acquire TCGA-STAD data from the GDC
#
# Downloads somatic mutations (MAF) and RNA-seq counts for stomach
# adenocarcinoma, then caches both as .rds so no later script has to hit the
# GDC API again. Downloads are slow; everything here is cache-guarded.
# ==============================================================================

source("00_setup.R")

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(maftools)
  library(dplyr)
})

PROJECT <- "TCGA-STAD"

# ---- Cohort overview ---------------------------------------------------------
gdc_projects <- getGDCprojects()
print(head(gdc_projects[, c("project_id", "name")]))
getProjectSummary(PROJECT)

# ---- Somatic mutations (masked MAF, open access) ------------------------------
maf <- cache_rds("stad_maf.rds", {
  snv_query <- GDCquery(
    project       = PROJECT,
    data.category = "Simple Nucleotide Variation",
    access        = "open",
    data.type     = "Masked Somatic Mutation",
    data.format   = "MAF"
  )
  GDCdownload(snv_query)
  read.maf(maf = GDCprepare(snv_query))
})

cat("Samples with mutation calls:", nrow(getSampleSummary(maf)), "\n")

# ---- RNA-seq counts (tumor + adjacent normal) --------------------------------
rna <- cache_rds("stad_rna.rds", {
  rna_query <- GDCquery(
    project               = PROJECT,
    data.category         = "Transcriptome Profiling",
    access                = "open",
    experimental.strategy = "RNA-Seq",
    sample.type           = c("Primary Tumor", "Solid Tissue Normal")
  )
  GDCdownload(rna_query)
  GDCprepare(rna_query)
})

cat("RNA-seq matrix:", nrow(rna), "genes x", ncol(rna), "samples\n")
print(table(colData(rna)$sample_type))

# ---- Gene annotation, split out for standalone reuse --------------------------
# Downstream scripts join on this to map Ensembl IDs to HGNC symbols.
gene_metadata <- cache_rds("gene_metadata.rds", as.data.frame(rowData(rna)))

# ---- Clinical / survival table -----------------------------------------------
clinical <- cache_rds("stad_clinical.rds",
                      GDCquery_clinic(project = PROJECT, type = "clinical"))

cat("Clinical records:", nrow(clinical), "\n")
