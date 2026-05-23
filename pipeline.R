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






























