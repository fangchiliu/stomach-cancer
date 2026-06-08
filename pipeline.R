# Project setup
getwd()


plots_dir = file.path(getwd(),"plots")
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

plots_dir = file.path(getwd(),"exports")
dir.create(plots_dir, 
           showWarnings = FALSE, 
           recursive = TRUE)

list.files()
list.dirs(getwd(), recursive = FALSE)

# Loading library
library(BiocManager)
library(remotes)
library(SummarizedExperiment)
library(TCGAbiolinks)
library(tidyverse)
library(maftools)
library(DESeq2)
library(EnhancedVolcano)
library(dplyr)
library(ggplot2)
library(ggrepel)

gdc_project = getGDCprojects()
head(gdc_project[c("project_id", "name")])
View(gdc_project)
gdc_project[c("project_id", "name")]
gdc_project$project_id
getProjectSummary("TCGA-STAD") 

snv_query = GDCquery(project = "TCGA-STAD",
                     data.category = "Simple Nucleotide Variation")

View(snv_query) #this is not readablefor us
output_snv_query=getResults(snv_query)
View(output_snv_query)

# unique(output_rna_query$sample_type)
# 
# GDCdownload(rna_query)
# 
# rna = GDCprepare(rna_query)

# Save the prepared RNA object
# saveRDS(rna, file = "exports/stad_rna.rds")
```


```{r}

# In a future/restarted session, load it back
rna <- readRDS("exports/stad_rna.rds")

# Extract expression matrix again
rna_matrix <- assay(rna, "unstranded")

View(rna_matrix)
dim(rna_matrix)
# head(rna_matrix)
```


differential expression analysis by using DESeq2
```{r}
# Create DESeq2 input object
dds_input <- DESeqDataSetFromMatrix(
  countData = rna_matrix,
  colData   = colData(rna),
  design    = ~ sample_type
)

# Check original factor levels
levels(dds_input$sample_type)

# Set Normal as reference level
# This makes log2FC = Tumor / Normal
dds_input$sample_type <- relevel(
  dds_input$sample_type,
  ref = "Solid Tissue Normal"
)

# Confirm new levels
levels(dds_input$sample_type)

# Filter low-count genes
keep <- rowSums(counts(dds_input)) >= 10
dds_input <- dds_input[keep, ]

# Run DESeq2 differential expression analysis
dds <- DESeq(dds_input)

# Extract results
res <- results(dds)

# View results
View(res)
res
summary(res)
res_df = as.data.frame(res)
View(res_df)


write.csv(res_df, file = 'exports/DEA.csv')

#Gene_id transfer
gene_metadata = as.data.frame(rowData(rna))
res_df$gene_id = rownames(res_df)
res_df = left_join(x = res_df,
                   y = gene_metadata,
                   by = 'gene_id')

dim(res_df)

EnhancedVolcano(toptable = res_df,
                lab = res_df$gene_name,
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1)



png(
  filename = file.path(plots_dir,"first_volacano.png"),
  width = 2600,
  height = 2000,
  res = 300
)
EnhancedVolcano(toptable = res_df,
                lab = res_df$gene_name,
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.01,
                FCcutoff = 3)
dev.off()

```



































