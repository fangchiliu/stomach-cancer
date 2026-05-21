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










