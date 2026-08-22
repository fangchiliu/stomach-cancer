# ==============================================================================
# 09_scrnaseq_analysis.R — Single-cell DE, QC, integration and UMAP
#
# Picks up the bundle from 08_scrnaseq_preprocessing.R and runs:
#   1. Compartment-wise cancer vs normal DE (presto Wilcoxon)
#   2. Seurat QC and filtering
#   3. Normalisation, HVGs, PCA
#   4. Restriction to Kumar et al. 2022 (the one study with both groups)
#   5. Harmony integration on donor, clustering, 2D + 3D UMAP
# ==============================================================================

source("00_setup.R")

suppressPackageStartupMessages({
  library(Matrix)
  library(tidyverse)
  library(Seurat)
  library(harmony)
  library(presto)      # remotes::install_github("immunogenomics/presto")
  library(hdf5r)
  library(patchwork)
  library(plotly)
  library(htmlwidgets)
})

H5AD_PATH <- file.path(base_dir, "scRNAseq_GC_N.h5ad")
N_DIMS    <- 30      # from the PCA elbow
STUDY     <- "Kumar et al. 2022"

bundle <- readRDS(exp_path("cancer_vs_normal_bundle.rds"))
mat    <- bundle$mat
meta   <- bundle$meta
stopifnot(all(colnames(mat) == meta$cell_id))

# ==============================================================================
# 1. Differential expression within each compartment
# ==============================================================================
run_de_compartment <- function(compartment_label) {
  cols <- which(meta$Compartment == compartment_label)
  grp  <- factor(meta$disease[cols], levels = c("normal", "gastric cancer"))

  cat("\n", compartment_label, "|", length(cols), "cells (",
      sum(grp == "normal"), "normal,", sum(grp == "gastric cancer"), "cancer )\n")

  # presto returns one row per (gene, group); the cancer rows are already
  # oriented as cancer vs the rest, which here means cancer vs normal.
  de <- as_tibble(wilcoxauc(mat[, cols, drop = FALSE], grp)) %>%
    filter(group == "gastric cancer") %>%
    transmute(
      gene          = feature,
      avgExpr_GC    = avgExpr,
      logFC_GC_vs_N = logFC,
      auc           = auc,          # > 0.5 means higher in cancer
      pct_expr_GC   = pct_in  / 100,
      pct_expr_N    = pct_out / 100,
      p_val         = pval,
      p_val_adj     = padj
    ) %>%
    arrange(p_val_adj, desc(abs(logFC_GC_vs_N)))

  tag <- gsub("[^A-Za-z]", "", compartment_label)
  write_csv(de, exp_path(sprintf("DE_%s_GC_vs_N_full.csv", tag)))

  # Require the gene to be detected somewhere, or tiny-percentage noise
  # dominates the top of the list.
  de_sig <- de %>%
    filter(!is.na(p_val_adj), p_val_adj < PADJ_CUTOFF,
           abs(logFC_GC_vs_N) > 0.25,
           pct_expr_GC > 0.05 | pct_expr_N > 0.05)

  write_csv(de_sig, exp_path(sprintf("DE_%s_GC_vs_N_sig.csv", tag)))
  cat("  significant genes:", nrow(de_sig), "\n")

  cat("\n  Top 10 up in cancer:\n")
  print(de %>% arrange(desc(logFC_GC_vs_N), p_val_adj) %>%
          dplyr::select(gene, logFC_GC_vs_N, auc, p_val_adj) %>% head(10))

  cat("\n  Top 10 down in cancer:\n")
  print(de %>% arrange(logFC_GC_vs_N, p_val_adj) %>%
          dplyr::select(gene, logFC_GC_vs_N, auc, p_val_adj) %>% head(10))

  invisible(de)
}

de_epithelial     <- run_de_compartment("Epithelial")
de_non_epithelial <- run_de_compartment("Non-Epithelial")

rm(mat, bundle); gc()

# ==============================================================================
# 2. Seurat object and QC
# ==============================================================================
counts  <- readRDS(exp_path("counts_cancer_vs_normal.rds"))
meta_df <- as.data.frame(meta)
rownames(meta_df) <- meta_df$cell_id
stopifnot(all(colnames(counts) == rownames(meta_df)))

so <- CreateSeuratObject(counts = counts, meta.data = meta_df,
                         project = "GC_vs_N", min.cells = 3, min.features = 0)

so[["percent.mt"]]   <- PercentageFeatureSet(so, pattern = "^MT-")
so[["percent.ribo"]] <- PercentageFeatureSet(so, pattern = "^RP[SL]")

print(so)
print(so@meta.data %>%
        group_by(disease) %>%
        summarise(n = n(),
                  med_nFeature = median(nFeature_RNA),
                  med_nCount   = median(nCount_RNA),
                  med_pct_mt   = round(median(percent.mt), 2),
                  .groups = "drop"))

p_vln <- VlnPlot(so, features = c("nFeature_RNA", "nCount_RNA",
                                  "percent.mt", "percent.ribo"),
                 group.by = "disease", pt.size = 0, ncol = 4)
ggsave(plot_path("QC_violin_by_disease.png"), p_vln,
       width = 16, height = 5, dpi = 150)

p_scatter <- FeatureScatter(so, "nCount_RNA", "percent.mt",   group.by = "disease") +
             FeatureScatter(so, "nCount_RNA", "nFeature_RNA", group.by = "disease")
ggsave(plot_path("QC_scatter.png"), p_scatter, width = 14, height = 6, dpi = 150)

# ==============================================================================
# 3. QC filtering
#
# Gentle thresholds: the atlas is already filtered upstream, so this mainly
# catches empty droplets and stressed cells.
# ==============================================================================
MIN_FEATURES <- 200
MAX_FEATURES <- 8000
MAX_MT       <- 20

md   <- so@meta.data
keep <- md$nFeature_RNA >= MIN_FEATURES &
        md$nFeature_RNA <= MAX_FEATURES &
        md$percent.mt   <  MAX_MT

cat("Cells:", ncol(so), "-> kept:", sum(keep),
    sprintf("(removed %.2f%%)\n", 100 * mean(!keep)))

# Check the filter isn't preferentially discarding one group or cell type,
# which would bias every downstream comparison.
cat("\nRemoval rate by disease:\n")
print(md %>% mutate(removed = !keep) %>% group_by(disease) %>%
        summarise(n = n(), pct_removed = round(100 * mean(removed), 2),
                  .groups = "drop"))

cat("\nHighest removal rates by cell type:\n")
print(md %>% mutate(removed = !keep) %>% group_by(Celltypes_global) %>%
        summarise(n = n(), pct_removed = round(100 * mean(removed), 2),
                  .groups = "drop") %>%
        arrange(desc(pct_removed)) %>% head(10))

so <- so[, keep]

# ==============================================================================
# 4. Normalise, HVGs, PCA, and attach batch metadata
# ==============================================================================
so <- NormalizeData(so, normalization.method = "LogNormalize", scale.factor = 1e4)
so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = 2000)

# Scale HVGs only — scaling all genes across 200k+ cells is not worth the memory.
so <- ScaleData(so, features = VariableFeatures(so))
so <- RunPCA(so, features = VariableFeatures(so), npcs = 50, verbose = FALSE)

ggsave(plot_path("PCA_elbow.png"), ElbowPlot(so, ndims = 50),
       width = 7, height = 5, dpi = 150)

# ---- Pull batch/donor columns back out of the h5ad ---------------------------------
h5  <- H5File$new(H5AD_PATH, mode = "r")
on.exit(try(h5$close_all(), silent = TRUE), add = TRUE)
obs <- h5[["obs"]]

read_h5ad_vector <- function(group, col) {
  obj <- group[[col]]
  if (inherits(obj, "H5Group")) {
    codes <- as.integer(obj[["codes"]]$read())
    cats  <- as.character(obj[["categories"]]$read())
    out <- rep(NA_character_, length(codes))
    ok  <- !is.na(codes) & codes >= 0
    out[ok] <- cats[codes[ok] + 1]
    return(out)
  }
  v <- as.character(obj$read()); v[v == ""] <- NA_character_; v
}

batch_cols <- intersect(
  c("Study", "Sample", "Batch", "donor_id",
    "Patient_status", "Tissue_in_paper", "tissue", "Location"),
  names(obs)
)

ridx <- so@meta.data$row_index_in_obs
for (cc in batch_cols) so[[cc]] <- read_h5ad_vector(obs, cc)[ridx]

cat("\nStudy x disease:\n")
print(table(so$Study, so$disease, useNA = "ifany"))

ggsave(
  plot_path("PCA_by_study_disease.png"),
  DimPlot(so, reduction = "pca", group.by = "Study",   raster = TRUE) + ggtitle("PCA by Study") +
  DimPlot(so, reduction = "pca", group.by = "disease", raster = TRUE) + ggtitle("PCA by disease"),
  width = 14, height = 6, dpi = 150
)

saveRDS(so, exp_path("so_GC_vs_N_pca.rds"))

# ==============================================================================
# 5. Restrict to a single study
#
# Across the full atlas, study and disease are confounded. Kumar et al. 2022 is
# the one study contributing both tumors and normal-adjacent tissue under the
# same protocol, so the contrast is run inside it.
# ==============================================================================
so_k <- subset(so, subset = Study == STUDY)
cat("\n", STUDY, "subset:", ncol(so_k), "cells\n")
print(table(so_k$disease))

donors_tumor  <- unique(so_k$donor_id[so_k$disease == "gastric cancer"])
donors_normal <- unique(so_k$donor_id[so_k$disease == "normal"])
cat("Tumor donors:", length(donors_tumor),
    "| normal donors:", length(donors_normal),
    "| paired:", length(intersect(donors_tumor, donors_normal)), "\n")

so_k <- FindVariableFeatures(so_k, selection.method = "vst", nfeatures = 2000)
so_k <- ScaleData(so_k, features = VariableFeatures(so_k))
so_k <- RunPCA(so_k, features = VariableFeatures(so_k), npcs = 50, verbose = FALSE)

ggsave(plot_path("Kumar_PCA_elbow.png"), ElbowPlot(so_k, ndims = 50),
       width = 7, height = 5, dpi = 150)

# ==============================================================================
# 6. Harmony integration, clustering, UMAP
# ==============================================================================
# Integrate on donor: patient identity is the dominant batch effect within a
# single study.
so_k <- RunHarmony(so_k, group.by.vars = "donor_id",
                   reduction.use = "pca", dims.use = 1:N_DIMS,
                   reduction.save = "harmony", verbose = FALSE)

so_k <- FindNeighbors(so_k, reduction = "harmony", dims = 1:N_DIMS, verbose = FALSE)
so_k <- FindClusters(so_k, resolution = 0.5, verbose = FALSE)
so_k <- RunUMAP(so_k, reduction = "harmony", dims = 1:N_DIMS, verbose = FALSE)

cat("\nClusters:", nlevels(so_k$seurat_clusters), "\n")
print(table(so_k$seurat_clusters, so_k$disease))

ggsave(
  plot_path("Kumar_UMAP_clusters_disease.png"),
  (DimPlot(so_k, group.by = "seurat_clusters", label = TRUE, raster = TRUE) +
     ggtitle("Clusters") + NoLegend()) +
  (DimPlot(so_k, group.by = "disease", raster = TRUE) + ggtitle("Disease")),
  width = 14, height = 6, dpi = 150
)

ggsave(
  plot_path("Kumar_UMAP_celltype_donor.png"),
  (DimPlot(so_k, group.by = "Celltypes_global", label = TRUE, repel = TRUE,
           raster = TRUE) + ggtitle("Cell type") + NoLegend()) +
  (DimPlot(so_k, group.by = "donor_id", raster = TRUE) + ggtitle("Donor") + NoLegend()),
  width = 14, height = 6, dpi = 150
)

saveRDS(so_k, exp_path("so_Kumar_clustered.rds"))

# ==============================================================================
# 7. Interactive 3D UMAP for the project site
# ==============================================================================
MAX_CELLS <- 45000   # three colour views triple the point count in the browser

if (!"umap3d" %in% Reductions(so_k)) {
  so_k <- RunUMAP(so_k, reduction = "harmony", dims = 1:N_DIMS, n.components = 3,
                  reduction.name = "umap3d", reduction.key = "UMAP3d_",
                  verbose = FALSE)
  saveRDS(so_k, exp_path("so_Kumar_clustered.rds"))
}

emb <- Embeddings(so_k, "umap3d")
df <- data.frame(
  UMAP1    = emb[, 1], UMAP2 = emb[, 2], UMAP3 = emb[, 3],
  CellType = as.character(so_k$Celltypes_global),
  disease  = as.character(so_k$disease),
  cluster  = as.character(so_k$seurat_clusters),
  stringsAsFactors = FALSE
)

# Downsample within cell type so rare populations survive.
if (nrow(df) > MAX_CELLS) {
  set.seed(1)
  df <- df %>% group_by(CellType) %>%
    slice_sample(prop = MAX_CELLS / nrow(df)) %>% ungroup()
}

# Maximally distinct categorical palette (Polychrome-style) — a gradient would
# make neighbouring cell types indistinguishable.
PALETTE <- c(
  "#AA0DFE","#3283FE","#85660D","#782AB6","#1C8356","#16FF32","#FE00FA",
  "#C4451C","#DEA0FD","#F7E1A0","#1CBE4F","#325A9B","#FEAF16","#F8A19F",
  "#90AD1C","#F6222E","#1CFFCE","#2ED9FF","#B10DA1","#C075A6","#FC1CBF",
  "#B00068","#FBE426","#FA0087","#00B5F7","#BDCDFF","#AAF400","#B5EFB5",
  "#7ED7D1","#66B0FF","#D85FF7","#FF6E3A","#FFD700","#00C2A0","#FF4FB6",
  "#8ADD00"
)

# Order cell types by compartment so the legend groups sensibly.
epi <- c("Foveolar_Differentiated","Foveolar_Intermediate","Neck-Cells","Chief",
         "Parietal","Enterocytes","Enterocytes_Intermediate",
         "Intestinal_Undifferentiated","Goblet","Enteroendocrine_CHGA",
         "Enteroendocrine_GAST","Enteroendocrine_GHRL","SMG","Superficial",
         "Suprabasal","Basal","Intermediate","Cancer")
imm <- c("T-Cells","Plasma-Cells","B-Cells","Macrophages","Mast")
str <- c("Fibroblasts","Endothelial")
bld <- c("Erythrocyte")

ct_order <- intersect(c(epi, imm, str, bld), unique(df$CellType))
ct_order <- c(ct_order, setdiff(unique(df$CellType), ct_order))
clusters <- as.character(sort(unique(as.integer(df$cluster))))

views <- list(
  CellType = list(var = "CellType", levels = ct_order, title = "Cell type",
                  pal = setNames(PALETTE[seq_along(ct_order)], ct_order)),
  Disease  = list(var = "disease", levels = c("gastric cancer", "normal"),
                  title = "Disease",
                  pal = c("gastric cancer" = "#E64B35", "normal" = "#4DBBD5")),
  Cluster  = list(var = "cluster", levels = clusters, title = "Cluster",
                  pal = setNames(rep(PALETTE, length.out = length(clusters)), clusters))
)
default_view <- "CellType"

# One trace per category per view; the dropdown toggles trace visibility.
p <- plot_ly()
trace_view <- character(0)

for (vn in names(views)) {
  v <- views[[vn]]
  for (lvl in v$levels) {
    sub <- df[df[[v$var]] == lvl, ]
    if (nrow(sub) == 0) next
    p <- add_trace(
      p, data = sub, x = ~UMAP1, y = ~UMAP2, z = ~UMAP3,
      type = "scatter3d", mode = "markers",
      name = lvl, legendgroup = vn,
      marker = list(size = 2.2, opacity = 0.7,
                    color = unname(v$pal[lvl]), line = list(width = 0)),
      visible = (vn == default_view), showlegend = TRUE,
      hovertemplate = paste0("<b>", lvl, "</b><extra></extra>")
    )
    trace_view <- c(trace_view, vn)
  }
}

bg   <- "#0e1117"
fg   <- "#e6e6e6"
grid <- "rgba(255,255,255,0.06)"
ax   <- list(title = "", showticklabels = FALSE, showgrid = TRUE,
             gridcolor = grid, zeroline = FALSE,
             showbackground = FALSE, showspikes = FALSE)

buttons <- lapply(names(views), function(vn) {
  list(method = "update", label = views[[vn]]$title,
       args = list(list(visible = trace_view == vn),
                   list("legend.title.text" = views[[vn]]$title)))
})

p <- p %>%
  layout(
    title = list(text = paste(STUDY, "— gastric 3D UMAP"),
                 font = list(color = fg, size = 18)),
    paper_bgcolor = bg, plot_bgcolor = bg,
    font = list(color = fg, family = "Inter, Helvetica, Arial, sans-serif"),
    legend = list(itemsizing = "constant", font = list(size = 11),
                  bgcolor = "rgba(0,0,0,0)", title = list(text = "Cell type")),
    scene = list(xaxis = ax, yaxis = ax, zaxis = ax, bgcolor = bg,
                 aspectmode = "data",
                 camera = list(eye = list(x = 1.5, y = 1.5, z = 0.9))),
    updatemenus = list(list(
      type = "dropdown", direction = "down", x = 0.01, y = 0.99,
      xanchor = "left", yanchor = "top", active = 0,
      bgcolor = "#1b2029", bordercolor = grid,
      font = list(color = fg, size = 12), buttons = buttons
    ))
  ) %>%
  config(responsive = TRUE, displaylogo = FALSE,
         toImageButtonOptions = list(format = "png", scale = 3,
                                     filename = "Kumar_UMAP3D"))

# Ship only the plotly.js modules actually used, for a much smaller HTML file.
p <- plotly::partial_bundle(p)
saveWidget(p, exp_path("Kumar_UMAP3D_colorby.html"),
           selfcontained = TRUE, title = "Kumar 3D UMAP")

# ---- Flat coordinate export for the website --------------------------------------
coords <- data.frame(
  cell_id     = colnames(so_k),
  UMAP1       = emb[, 1], UMAP2 = emb[, 2], UMAP3 = emb[, 3],
  CellType    = as.character(so_k$Celltypes_global),
  disease     = as.character(so_k$disease),
  cluster     = as.character(so_k$seurat_clusters),
  compartment = as.character(so_k$Compartment_fine),
  donor       = as.character(so_k$donor_id),
  stringsAsFactors = FALSE
)

if ("umap" %in% Reductions(so_k)) {
  emb2d <- Embeddings(so_k, "umap")
  coords$UMAP2D_1 <- emb2d[, 1]
  coords$UMAP2D_2 <- emb2d[, 2]
}

write_csv(coords, exp_path("Kumar_UMAP_coords_labels.csv"))
write_csv(data.frame(CellType = ct_order, color = PALETTE[seq_along(ct_order)]),
          exp_path("Kumar_celltype_colors.csv"))

cat("Single-cell analysis complete.\n")
