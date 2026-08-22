# ==============================================================================
# 08_scrnaseq_preprocessing.R — Extract gastric cancer / normal cells from h5ad
#
# The CellxGene single-cell atlas ships as one large .h5ad. Rather than load it
# through anndata, this reads the sparse matrix directly with hdf5r and pulls
# only the cells we need, which keeps memory manageable.
#
# Output: a log-normalised gene x cell matrix plus tagged metadata, bundled
# into one .rds that 09_scrnaseq_analysis.R reloads.
# ==============================================================================

source("00_setup.R")

suppressPackageStartupMessages({
  library(hdf5r)
  library(Matrix)
  library(tidyverse)
})

H5AD_PATH   <- file.path(base_dir, "scRNAseq_GC_N.h5ad")
KEEP_STATES <- c("normal", "gastric cancer")

stopifnot(file.exists(H5AD_PATH))

# ==============================================================================
# Helpers for reading AnnData's HDF5 layout
# ==============================================================================

# Categorical obs/var columns are stored as integer codes plus a category
# table; older AnnData versions kept categories in a sibling __categories group.
read_h5ad_vector <- function(group, col) {
  obj <- group[[col]]

  decode <- function(codes, cats) {
    out <- rep(NA_character_, length(codes))
    ok  <- !is.na(codes) & codes >= 0
    out[ok] <- cats[codes[ok] + 1]
    out
  }

  if (inherits(obj, "H5Group")) {
    codes <- as.integer(obj[["codes"]]$read())
    if ("categories" %in% names(obj)) {
      return(decode(codes, as.character(obj[["categories"]]$read())))
    }
    if ("__categories" %in% names(group) && col %in% names(group[["__categories"]])) {
      return(decode(codes, as.character(group[["__categories"]][[col]]$read())))
    }
    stop("Unknown encoding for column: ", col)
  }

  v <- as.character(obj$read())
  v[v == ""] <- NA_character_
  v
}

# Read an arbitrary subset of rows out of a CSR matrix without materialising
# the whole thing as dense.
read_h5ad_csr_rows <- function(h5, x_path, row_idx, n_cols) {
  xg <- h5
  for (p in strsplit(x_path, "/", fixed = TRUE)[[1]]) xg <- xg[[p]]

  data    <- xg[["data"]]$read()
  indices <- as.integer(xg[["indices"]]$read())   # 0-based column indices
  indptr  <- as.integer(xg[["indptr"]]$read())    # length n_rows + 1

  row_idx_0 <- as.integer(row_idx - 1L)
  stopifnot(all(row_idx_0 >= 0), all(row_idx_0 < length(indptr) - 1L))

  i_list <- j_list <- x_list <- vector("list", length(row_idx_0))

  for (k in seq_along(row_idx_0)) {
    r     <- row_idx_0[k]
    start <- indptr[r + 1L] + 1L
    end   <- indptr[r + 2L]

    if (end < start) {
      i_list[[k]] <- integer(0); j_list[[k]] <- integer(0); x_list[[k]] <- numeric(0)
      next
    }

    cols0 <- indices[start:end]
    i_list[[k]] <- rep.int(k, length(cols0))
    j_list[[k]] <- as.integer(cols0 + 1L)
    x_list[[k]] <- as.numeric(data[start:end])
  }

  sparseMatrix(i = unlist(i_list), j = unlist(j_list), x = unlist(x_list),
               dims = c(length(row_idx), n_cols), repr = "R")
}

# Sum rows sharing a name, staying sparse throughout.
collapse_duplicate_rows <- function(mat) {
  if (!any(duplicated(rownames(mat)))) return(mat)
  fac <- factor(rownames(mat), levels = unique(rownames(mat)))
  agg <- sparseMatrix(i = as.integer(fac), j = seq_along(fac), x = 1,
                      dims = c(nlevels(fac), length(fac)))
  out <- agg %*% mat
  rownames(out) <- levels(fac)
  as(out, "dgCMatrix")
}

# ==============================================================================
# 1. Inspect obs, then subset to the two disease states
# ==============================================================================
h5 <- H5File$new(H5AD_PATH, mode = "r")
on.exit(try(h5$close_all(), silent = TRUE), add = TRUE)

obs <- h5[["obs"]]
var <- h5[["var"]]

idx_name <- if ("_index" %in% names(obs)) "_index" else "__index__"
cell_ids <- as.character(obs[[idx_name]]$read())

disease            <- read_h5ad_vector(obs, "disease")
celltypes_global   <- read_h5ad_vector(obs, "Celltypes_global")
detailed_cell_type <- read_h5ad_vector(obs, "Detailed_Cell_Type")

keep <- !is.na(disease) & disease %in% KEEP_STATES

meta <- tibble(
  cell_id            = cell_ids[keep],
  disease            = disease[keep],
  Celltypes_global   = celltypes_global[keep],
  Detailed_Cell_Type = detailed_cell_type[keep],
  row_index_in_obs   = which(keep)
)

cat("Cells kept (cancer + normal):", nrow(meta), "\n")
print(table(meta$disease))
print(table(meta$Celltypes_global, meta$disease))

# ==============================================================================
# 2. Tag each cell with a tissue compartment
#
# The DE contrast is run within compartment, so a tumor-vs-normal difference
# can't just be a shift in cell type proportions.
# ==============================================================================
epithelial_types <- c(
  "Foveolar_Differentiated", "Foveolar_Intermediate", "Neck-Cells",
  "Chief", "Parietal", "Enterocytes", "Enterocytes_Intermediate",
  "Intestinal_Undifferentiated", "Goblet",
  "Enteroendocrine_CHGA", "Enteroendocrine_GAST", "Enteroendocrine_GHRL",
  "SMG", "Superficial", "Suprabasal", "Basal", "Intermediate", "Cancer"
)
immune_types <- c("T-Cells", "Plasma-Cells", "Macrophages", "B-Cells", "Mast")
stroma_types <- c("Fibroblasts", "Endothelial")
blood_types  <- c("Erythrocyte")

meta <- meta %>%
  mutate(
    Compartment_fine = case_when(
      Celltypes_global %in% epithelial_types ~ "Epithelial",
      Celltypes_global %in% immune_types     ~ "Immune",
      Celltypes_global %in% stroma_types     ~ "Stroma",
      Celltypes_global %in% blood_types      ~ "Blood",
      TRUE                                   ~ "UNCLASSIFIED"
    ),
    Compartment = case_when(
      Compartment_fine == "Epithelial"                     ~ "Epithelial",
      Compartment_fine %in% c("Immune", "Stroma", "Blood") ~ "Non-Epithelial",
      TRUE                                                 ~ "UNCLASSIFIED"
    )
  )

unclassified <- meta %>% filter(Compartment == "UNCLASSIFIED") %>% count(Celltypes_global)
if (nrow(unclassified) > 0) {
  warning("Unclassified cell types found:")
  print(unclassified)
}

print(table(meta$Compartment_fine, meta$disease))
write_csv(meta, exp_path("meta_cancer_vs_normal_tagged.csv"))

# ==============================================================================
# 3. Read the expression matrix for those cells
# ==============================================================================
var_idx_name <- if ("_index" %in% names(var)) "_index" else "__index__"
gene_ids     <- as.character(var[[var_idx_name]]$read())

# raw/X holds counts when the atlas ships pre-normalised X; prefer it so we
# control the normalisation ourselves.
x_path <- if ("raw" %in% names(h5) && "X" %in% names(h5[["raw"]])) "raw/X" else "X"
cat("Reading matrix from:", x_path, "|", nrow(meta), "cells\n")

counts_cg <- read_h5ad_csr_rows(h5, x_path, meta$row_index_in_obs, length(gene_ids))
rownames(counts_cg) <- meta$cell_id
colnames(counts_cg) <- gene_ids

counts <- as(as(t(counts_cg), "CsparseMatrix"), "dgCMatrix")   # genes x cells

# Ensembl IDs -> symbols; duplicated symbols get their counts summed.
feature_name <- read_h5ad_vector(var, "feature_name")
rownames(counts) <- ifelse(is.na(feature_name) | feature_name == "",
                           rownames(counts), feature_name)
counts <- collapse_duplicate_rows(counts)

cat("Counts matrix:", nrow(counts), "genes x", ncol(counts), "cells\n")

# ---- Log-normalise ----------------------------------------------------------------
mat <- counts
if (x_path == "raw/X") {
  sf     <- 1e4 / pmax(Matrix::colSums(mat), 1)
  mat    <- mat %*% Diagonal(x = sf)
  mat@x  <- log1p(mat@x)
  mat    <- as(mat, "dgCMatrix")
  cat("Applied log1p(CP10K) normalization.\n")
}

# ---- Marker sanity check ----------------------------------------------------------
markers <- c("EPCAM", "KRT8", "KRT18", "KRT19", "MUC5AC", "MUC6",
             "TFF1", "TFF2", "TFF3", "PTPRC", "CD3D", "CD68",
             "COL1A1", "PECAM1", "CDX2")
print(data.frame(gene = markers, present = markers %in% rownames(mat)))

# ==============================================================================
# 4. Save aligned bundle
# ==============================================================================
stopifnot(all(colnames(mat) == meta$cell_id))

saveRDS(counts, exp_path("counts_cancer_vs_normal.rds"))

bundle <- list(
  mat  = mat,
  meta = meta,
  info = list(
    source_h5ad      = basename(H5AD_PATH),
    diseases         = KEEP_STATES,
    normalization    = "log1p(CP10K) from raw/X",
    dims_genes_cells = dim(mat),
    created          = Sys.time()
  )
)

saveRDS(bundle, exp_path("cancer_vs_normal_bundle.rds"))

cat("Saved bundle:", dim(mat)[1], "genes x", dim(mat)[2], "cells |",
    format(object.size(bundle), units = "auto"), "\n")
