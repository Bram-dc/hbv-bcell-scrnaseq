# Re-cluster an antigen subset (Ig genes dropped from HVG, not from the object); UMAP is computed last so it cannot perturb cluster assignments.

suppressPackageStartupMessages(library(Seurat))

cluster_umap <- function(obj, dims = PCA_DIMS, res = CLUSTER_RES,
                         drop_ig = TRUE, with_umap = TRUE) {
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, verbose = FALSE)
  if (drop_ig) {
    ig <- intersect(ig_genes(), rownames(obj))
    VariableFeatures(obj) <- setdiff(VariableFeatures(obj), ig)
  }
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = dims, verbose = FALSE)
  obj <- FindClusters(obj, resolution = res, verbose = FALSE)
  if (with_umap) {
    obj <- RunUMAP(obj, dims = dims, verbose = FALSE)
  }
  obj
}

label_clusters <- function(obj, antigen = NULL) {
  cl <- as.character(obj$seurat_clusters)
  lev <- as.character(sort(unique(as.integer(cl))))
  labs <- if (!is.null(antigen)) CLUSTER_LABELS[[antigen]] else NULL
  named <- vapply(lev, function(k) {
    if (!is.null(labs) && k %in% names(labs)) labs[[k]] else k
  }, character(1))
  obj$cluster_label <- factor(unname(named[cl]), levels = unique(unname(named)))
  obj
}
