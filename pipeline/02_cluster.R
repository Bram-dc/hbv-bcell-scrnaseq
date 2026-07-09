# Per-antigen clustering EDA, a clustree resolution sweep, and the IA/IC/FC HBcAg embedding overlay. -> output/clusters/<antigen>/

source("config/config.R")
source("config/comparisons.R")
invisible(lapply(list.files("R", "[.]R$", full.names = TRUE), source))

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(clustree)
})

ANTIGENS <- c("HBcAg")
RESOLUTIONS <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2)

obj <- prepare()

for (antigen in ANTIGENS) {
  message(sprintf("\n== Cluster EDA: %s ==", antigen))
  out_dir <- file.path(OUTPUT_ROOT, "clusters", antigen)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  ag <- label_clusters(cluster_umap(subset_antigen(obj, antigen)), antigen)

  umap_cluster_plot(ag, file.path(out_dir, "umap_clusters.png"),
    title = sprintf("%s-specific B cells - unsupervised clusters", antigen)
  )
  umap_meta_plot(ag, BATCH_VAR, file.path(out_dir, "umap_sequencing_batch.png"),
    title = sprintf("%s-specific B cells - sequencing batch", antigen)
  )
  composition_plot(ag, file.path(out_dir, "cluster_composition_per_group.png"))
  if (antigen == "HBcAg") {
    composition_state_plot(ag, file.path(out_dir, "cluster_composition_per_state.png"))
    cohort_state_plot(ag, file.path(out_dir, "cohort_patients_cells_by_state.png"))
  }
  qc_violin_plot(ag, file.path(out_dir, "qc_metrics_per_cluster.png"))
  marker_dotplot(ag, ALL_MARKERS, file.path(out_dir, "dotplot_markers_per_cluster.png"))
  violin_markers_plot(ag, VIOLIN_MARKERS, file.path(out_dir, "violin_markers_per_cluster.png"))
  cluster_markers(ag, out_dir)

  # Clustree resolution sweep: re-cluster the antigen subset across resolutions.
  cl <- NormalizeData(subset_antigen(obj, antigen), verbose = FALSE)
  cl <- FindVariableFeatures(cl, verbose = FALSE)
  VariableFeatures(cl) <- setdiff(VariableFeatures(cl), intersect(ig_genes(), rownames(cl)))
  cl <- ScaleData(cl, verbose = FALSE)
  cl <- RunPCA(cl, verbose = FALSE)
  cl <- FindNeighbors(cl, dims = PCA_DIMS, verbose = FALSE)
  for (res in RESOLUTIONS) cl <- FindClusters(cl, resolution = res, verbose = FALSE)
  p_tree <- clustree(cl, prefix = "RNA_snn_res.") +
    labs(title = "Cluster stability across resolutions", subtitle = antigen)
  ggsave(file.path(out_dir, "clustree.png"), p_tree,
    width = 12, height = 14,
    dpi = 200, bg = "white"
  )
}

# Shared HBcAg embedding is valid because clusters are stable across contrasts.
hb <- cluster_umap(subset_antigen(obj, "HBcAg"))
gc <- trimws(hb$Group)
hb$condition <- dplyr::case_when(
  gc == "IA" ~ "IA",
  gc %in% c("IC(high)", "IC(low)") ~ "IC",
  gc %in% FC_GROUPS ~ "FC", TRUE ~ "Other"
)
ud <- as.data.frame(Embeddings(hb, "umap"))
colnames(ud) <- c("UMAP_1", "UMAP_2")
ud$Cluster <- factor(as.character(hb$seurat_clusters),
  levels = as.character(sort(as.integer(unique(as.character(hb$seurat_clusters)))))
)
ud$Condition <- factor(hb$condition, levels = c("IA", "IC", "FC", "Other"))
cond_pal <- c("IA" = "#E63946", "IC" = "#F4A261", "FC" = "#457B9D", "Other" = "grey80")
centroids <- ud |>
  group_by(Cluster) |>
  summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = "drop")
p_overlay <- ggplot() +
  geom_point(
    data = filter(ud, Condition == "Other"), aes(UMAP_1, UMAP_2),
    colour = "grey85", size = 0.8, alpha = 0.5
  ) +
  geom_point(
    data = filter(ud, Condition != "Other"), aes(UMAP_1, UMAP_2, colour = Condition),
    size = 1.2, alpha = 0.75
  ) +
  geom_text(data = centroids, aes(UMAP_1, UMAP_2, label = Cluster), fontface = "bold", size = 5) +
  scale_colour_manual(values = cond_pal, breaks = c("IA", "IC", "FC")) +
  labs(
    title = "IA, IC and FC on the shared HBcAg embedding",
    subtitle = "Grey = other groups (clustered together, not compared)"
  ) +
  theme_minimal(base_size = 12)
ggsave(file.path(OUTPUT_ROOT, "clusters", "HBcAg", "umap_condition_overlay.png"),
  p_overlay,
  width = 8, height = 6, dpi = 200, bg = "white"
)

message("\nCluster EDA complete.")
