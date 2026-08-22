# Gastric Cancer Multi-Omics Analysis

Identification and validation of transcriptional targets in stomach
adenocarcinoma, across bulk RNA-seq, microarray, and single-cell data.

The analysis converges on **INHBA**: significantly upregulated in tumors in two
independent cohorts, associated with worse overall survival, and co-regulated
with the EMT program.

## Data sources

| Cohort | Platform | Role |
|---|---|---|
| TCGA-STAD | RNA-seq + somatic MAF (GDC) | Discovery |
| GSE66229 (ACRG) | Affymetrix microarray (GEO) | Independent validation |
| CellxGene gastric atlas | 10x single-cell RNA-seq | Cell-type resolution |

## Pipeline

Scripts are numbered in dependency order. Each one sources `00_setup.R` and
caches its expensive steps to `exports/`, so a rerun skips work already done.

| Script | What it does |
|---|---|
| `00_setup.R` | Directory layout, shared constants and palette, helpers |
| `01_tcga_download.R` | Pulls TCGA-STAD mutations, RNA-seq counts and clinical data from the GDC |
| `02_mutation_landscape.R` | Cohort mutation burden, oncoplot, Ti/Tv spectrum |
| `03_differential_expression.R` | DESeq2 tumor vs normal; exports the VST matrix and volcano plot |
| `04_pathway_enrichment.R` | GSEA against Hallmark, Reactome, GO:BP and KEGG; EMT figure |
| `05_survival_analysis.R` | Kaplan-Meier overall survival by INHBA median split |
| `06_validation_gse66229.R` | limma DE in the ACRG cohort, then cross-cohort intersection |
| `07_inhba_characterization.R` | INHBA expression in both cohorts; co-expression-driven GSEA |
| `08_scrnaseq_preprocessing.R` | Reads the `.h5ad` atlas directly, subsets and tags cells by compartment |
| `09_scrnaseq_analysis.R` | Compartment-wise DE, QC, Harmony integration, 2D/3D UMAP |

## Running it

```r
setwd("path/to/project")

source("00_setup.R")   # creates exports/ and plots/

source("01_tcga_download.R")
source("02_mutation_landscape.R")
source("03_differential_expression.R")
source("04_pathway_enrichment.R")
source("05_survival_analysis.R")
source("06_validation_gse66229.R")
source("07_inhba_characterization.R")
```

The single-cell scripts need `scRNAseq_GC_N.h5ad` in the working directory:

```r
source("08_scrnaseq_preprocessing.R")
source("09_scrnaseq_analysis.R")
```

Package installation is listed (commented out) at the bottom of `00_setup.R`.

## Notes on method choices

- **Reference level.** `sample_type` is releveled to `Solid Tissue Normal`, so
  every log2 fold change reads as tumor over normal.
- **GSEA ranking metric.** Genes are ranked by
  `sign(log2FC) * -log10(p)` rather than fold change alone, which keeps both the
  direction and the strength of evidence.
- **Cross-cohort threshold.** The intersection uses `|log2FC| > 0.5` on both
  sides. Microarray intensities compress fold changes relative to RNA-seq, so a
  looser common cut keeps the comparison symmetric.
- **Survival censoring.** Late TCGA follow-up is sparse, so the headline
  Kaplan-Meier curve is administratively censored at 5 years.
- **Single-cell study restriction.** Across the full atlas, study and disease
  state are confounded. The tumor vs normal contrast is therefore run inside
  Kumar et al. 2022, the one study contributing both under a shared protocol.
- **Single-cell batch correction.** Harmony is run on `donor_id`, since patient
  identity is the dominant batch effect within a single study.

## Output layout

```
exports/    result tables (.csv), cached objects (.rds), interactive UMAP (.html)
plots/      figures (.png / .pdf)
```

Neither directory is tracked in git; both are recreated by `00_setup.R`.
