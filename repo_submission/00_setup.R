# ==============================================================================
# 00_setup.R — Shared project configuration
#
# Sourced by every downstream script. Defines directory layout, plotting
# constants, and small helpers used across the analysis.
# ==============================================================================

# ---- Directory layout --------------------------------------------------------
base_dir    <- getwd()
exports_dir <- file.path(base_dir, "exports")   # tables and cached R objects
plots_dir   <- file.path(base_dir, "plots")     # figures

dir.create(exports_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir,   showWarnings = FALSE, recursive = TRUE)

# ---- Analysis constants ------------------------------------------------------
PADJ_CUTOFF <- 0.05
LFC_CUTOFF  <- 1      # TCGA RNA-seq (log2)
LFC_ARRAY   <- 0.5    # microarray log2 ratios compress, so a looser cut
SEED        <- 42

# Shared palette so every figure in the project reads as one set.
COL_TUMOR   <- "#E63946"
COL_NORMAL  <- "#457B9D"
COL_NEUTRAL <- "grey80"
COL_ACCENT  <- "#F1A208"

# ---- Helpers -----------------------------------------------------------------

# Path builders, so scripts never hard-code "exports/" or "plots/".
exp_path  <- function(...) file.path(exports_dir, ...)
plot_path <- function(...) file.path(plots_dir, ...)

# Save a ggplot at print quality with a white background.
save_plot <- function(plot, filename, width = 9, height = 7, dpi = 300) {
  ggplot2::ggsave(plot_path(filename), plot,
                  width = width, height = height, dpi = dpi, bg = "white")
  message("Saved: ", plot_path(filename))
}

# Compute an expensive object once, then reuse the cached copy on later runs.
cache_rds <- function(filename, expr) {
  path <- exp_path(filename)
  if (file.exists(path)) {
    message("Loading cached: ", path)
    return(readRDS(path))
  }
  obj <- force(expr)
  saveRDS(obj, path)
  message("Cached: ", path)
  obj
}

# ---- Package installation (run once, interactively) --------------------------
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c(
#   "SummarizedExperiment", "TCGAbiolinks", "maftools", "DESeq2",
#   "EnhancedVolcano", "clusterProfiler", "org.Hs.eg.db", "enrichplot",
#   "pathview", "GEOquery", "limma", "Biobase"
# ))
# install.packages(c("tidyverse", "ggrepel", "ggpubr", "survival", "survminer",
#                    "msigdbr", "patchwork", "Seurat", "harmony", "hdf5r",
#                    "plotly", "htmlwidgets"))
# remotes::install_github("immunogenomics/presto")

theme_set_project <- function() {
  ggplot2::theme_set(
    ggplot2::theme_minimal(base_size = 13) +
      ggplot2::theme(
        plot.title       = ggplot2::element_text(face = "bold"),
        plot.subtitle    = ggplot2::element_text(color = "grey40"),
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_line(color = "grey93", linewidth = 0.3)
      )
  )
}
