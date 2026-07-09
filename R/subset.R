# Subsetting: antigen subsets, comparison resolution, and single-cluster slicing.

suppressPackageStartupMessages(library(Seurat))

# Cells of one or more antigens.
subset_antigen <- function(obj, antigen) {
  meta <- obj@meta.data
  cells <- rownames(meta)[trimws(meta$Antigen) %in% antigen]
  subset(obj, cells = cells)
}

# Cluster the full antigen subset before restricting to two groups, so cluster ids stay stable across contrasts.
subset_comparison <- function(obj, spec, dims = PCA_DIMS, res = CLUSTER_RES) {
  ag <- cluster_umap(subset_antigen(obj, spec$antigen), dims = dims, res = res)

  gc <- trimws(ag@meta.data[[spec$group_col]])
  cond <- ifelse(gc %in% spec$groups1, spec$ident1,
    ifelse(gc %in% spec$groups2, spec$ident2, NA_character_)
  )

  excl <- if (!is.null(spec$exclude_clusters)) spec$exclude_clusters else character(0)
  present <- unique(as.character(ag$seurat_clusters))
  missing <- setdiff(excl, present)
  if (length(missing) > 0) {
    warning(
      sprintf(
        "[%s] exclude_clusters %s not present in the %s clustering (clusters: %s) - ignored",
        spec$name, paste(missing, collapse = ", "),
        paste(spec$antigen, collapse = "+"), paste(sort(present), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  keep <- !is.na(cond) & !(as.character(ag$seurat_clusters) %in% excl)

  sub <- subset(ag, cells = colnames(ag)[keep])
  sub$condition <- factor(cond[keep], levels = c(spec$ident1, spec$ident2))
  sub
}

# Single-cluster slice.
subset_cluster <- function(obj, cl) {
  subset(obj, cells = colnames(obj)[as.character(obj$seurat_clusters) == cl])
}

# Alias used by the SCENIC driver (Rscript ... <comparison> cluster_0).
with_clusters <- function(obj, cluster_id) {
  cl <- sub("^cluster_", "", cluster_id)
  subset_cluster(obj, cl)
}

# Sorted character cluster ids of an object.
cluster_ids <- function(obj) {
  as.character(sort(as.integer(as.character(unique(obj$seurat_clusters)))))
}
