# Data loading and preparation into a QC'd, doublet-filtered, clustered Seurat object.

suppressPackageStartupMessages({
  library(Seurat)
  library(DoubletFinder)
})

# Load the raw object and clean the key metadata columns.
load_dataset <- function(rdata_path = DATA_FILE, object_name = OBJECT_NAME) {
  load(rdata_path)
  if (!exists(object_name)) {
    stop("Object '", object_name, "' not found in ", rdata_path)
  }
  obj <- get(object_name)
  obj$Group <- trimws(obj$Group)
  obj$Antigen <- trimws(obj$Antigen)
  obj
}

# Cell QC: nFeature 200-3000, percent.mt < 20, nCount > 500.
apply_qc <- function(obj) {
  n_before <- ncol(obj)
  obj <- subset(
    obj,
    subset = nFeature_RNA > QC_MIN_FEATURES &
      nFeature_RNA < QC_MAX_FEATURES &
      percent.mt < QC_MAX_MT &
      nCount_RNA > QC_MIN_COUNTS
  )
  message("QC filter: kept ", ncol(obj), " of ", n_before, " cells")
  obj
}

# Gene filter: keep protein-coding (org.Hs.eg.db) plus Ig (HGNC).
# Ig genes are often not typed "protein-coding", so the union retains them.
filter_genes <- function(obj) {
  n_before <- nrow(obj)
  coding <- protein_coding_genes(rownames(obj))
  ig <- tryCatch(ig_genes(), error = function(e) {
    warning("HGNC not available; Ig genes not added to keep-list")
    character(0)
  })
  keep <- intersect(rownames(obj), union(coding, ig))
  obj <- subset(obj, features = keep)
  message(
    "Gene filter: kept ", nrow(obj), " of ", n_before,
    " genes (protein-coding + Ig)"
  )
  obj
}

# DoubletFinder: keep singlets only.
remove_doublets <- function(obj, dims = PCA_DIMS) {
  annotations <- obj$seurat_clusters
  homotypic <- modelHomotypic(annotations)
  nExp_poi <- round(DF_NEXP_FRAC * ncol(obj))
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic))

  obj <- doubletFinder(
    obj,
    PCs = dims, pN = DF_PN, pK = DF_PK, nExp = nExp_poi.adj
  )
  singlet_col <- grep("^DF.classifications_", colnames(obj@meta.data), value = TRUE)
  if (length(singlet_col) == 0) {
    stop("No DoubletFinder classification column found in metadata")
  }
  singlet_col <- singlet_col[1]
  singlets <- rownames(obj@meta.data)[obj@meta.data[[singlet_col]] == "Singlet"]
  subset(obj, cells = singlets)
}

# Full preparation: load -> QC/gene filter -> normalize/PCA/cluster -> DoubletFinder -> cache.
# apply_qc = FALSE keeps all cells/genes for the QC-transition diagnostics (cached separately).
prepare <- function(apply_qc = TRUE,
                    use_cache = TRUE,
                    rdata_path = DATA_FILE,
                    object_name = OBJECT_NAME,
                    cache_path = if (apply_qc) SINGLET_CACHE else PREQC_CACHE) {
  if (use_cache && file.exists(cache_path)) {
    message("Loading cached Seurat object: ", cache_path)
    return(readRDS(cache_path))
  }

  obj <- load_dataset(rdata_path, object_name)

  if (apply_qc) {
    obj <- apply_qc_filter_genes(obj)
  } else {
    message("QC filter SKIPPED (apply_qc = FALSE): ", ncol(obj), " cells")
  }

  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj)
  obj <- FindNeighbors(obj, dims = PCA_DIMS)
  obj <- FindClusters(obj)
  obj <- remove_doublets(obj, dims = PCA_DIMS)

  dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(obj, cache_path)
  message("Saved cache: ", cache_path)
  obj
}

# Cell QC + gene filter together.
apply_qc_filter_genes <- function(obj) {
  filter_genes(apply_qc(obj))
}
