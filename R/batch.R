# Batch-effect diagnostics on the whole dataset: UMAP per batch, cluster composition, and cluster/PC-batch association.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

# Cramer's V for two categoricals: 0 = independent, 1 = fully determined.
cramers_v <- function(cat1, cat2) {
  tbl <- table(cat1, cat2)
  chi <- suppressWarnings(chisq.test(tbl))
  n <- sum(tbl)
  k <- min(nrow(tbl), ncol(tbl)) - 1
  v <- if (k > 0) sqrt((unname(chi$statistic) / n) / k) else NA_real_
  data.frame(cramers_v = v, p_value = chi$p.value)
}

batch_diagnostics <- function(obj, out_dir,
                              batch_candidates = c(
                                "Sequencing_batch",
                                "Processing_batch", "Sending_batch", "Library"
                              )) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  batch_vars <- intersect(batch_candidates, colnames(obj@meta.data))
  for (bv in batch_vars) {
    obj@meta.data[[bv]] <- trimws(as.character(obj@meta.data[[bv]]))
  }
  if (!"umap" %in% Reductions(obj)) {
    obj <- RunUMAP(obj, dims = PCA_DIMS, verbose = FALSE)
  }
  clusters <- obj$seurat_clusters

  umap_df <- as.data.frame(Embeddings(obj, "umap"))
  colnames(umap_df)[1:2] <- c("UMAP_1", "UMAP_2")
  umap_df <- cbind(umap_df, obj@meta.data[, batch_vars, drop = FALSE])

  batch_umaps <- lapply(batch_vars, function(bv) {
    ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, colour = .data[[bv]])) +
      geom_point(size = 0.3, alpha = 0.6) +
      guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
      labs(title = paste("UMAP coloured by", bv), colour = bv) +
      theme_minimal(base_size = 11) +
      theme(legend.key.size = unit(0.35, "cm"))
  })
  ggsave(file.path(out_dir, "batch_umap.png"), wrap_plots(batch_umaps, ncol = 2),
    width = 13, height = 10, dpi = 200, bg = "white"
  )

  comp_df <- bind_rows(lapply(batch_vars, function(bv) {
    obj@meta.data |>
      transmute(Cluster = clusters, Batch = .data[[bv]]) |>
      count(Cluster, Batch) |>
      group_by(Cluster) |>
      mutate(prop = n / sum(n)) |>
      ungroup() |>
      mutate(BatchVar = bv)
  }))
  p_comp <- ggplot(comp_df, aes(x = Cluster, y = prop, fill = Batch)) +
    geom_col(position = "stack", width = 0.85) +
    facet_wrap(~BatchVar, ncol = 1, scales = "free_y") +
    scale_y_continuous(labels = scales::percent) +
    labs(
      title = "Batch composition per cluster",
      subtitle = "Even mixing indicates low batch effect",
      x = "Seurat cluster", y = "Proportion"
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.key.size = unit(0.35, "cm"))
  ggsave(file.path(out_dir, "batch_cluster_comp.png"), p_comp,
    width = 11, height = 4 * length(batch_vars), dpi = 200, bg = "white"
  )

  assoc_df <- bind_rows(lapply(batch_vars, function(bv) {
    res <- cramers_v(clusters, obj@meta.data[[bv]])
    res$BatchVar <- bv
    res
  }))
  write.csv(assoc_df, file.path(out_dir, "batch_cluster_association.csv"), row.names = FALSE)
  p_assoc <- ggplot(assoc_df, aes(x = reorder(BatchVar, cramers_v), y = cramers_v, fill = cramers_v)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = sprintf("V=%.2f\np=%.1e", cramers_v, p_value)), hjust = -0.1, size = 3) +
    scale_fill_gradient(low = "#A8DADC", high = "#E63946", limits = c(0, 1)) +
    coord_flip() +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.25))) +
    labs(
      title = "Cluster vs batch association (Cramer's V)",
      subtitle = "Higher = clusters more determined by batch", x = NULL, y = "Cramer's V"
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none")
  ggsave(file.path(out_dir, "batch_cluster_assoc.png"), p_assoc,
    width = 8, height = 5, dpi = 200, bg = "white"
  )

  n_pc <- min(20, ncol(Embeddings(obj, "pca")))
  pca_emb <- Embeddings(obj, "pca")[, 1:n_pc, drop = FALSE]
  pc_assoc <- bind_rows(lapply(batch_vars, function(bv) {
    grp <- obj@meta.data[[bv]]
    data.frame(
      PC = factor(paste0("PC", 1:n_pc), levels = paste0("PC", 1:n_pc)),
      BatchVar = bv,
      neglog10p = vapply(1:n_pc, function(i) {
        p <- suppressWarnings(kruskal.test(pca_emb[, i], factor(grp))$p.value)
        -log10(pmax(p, 1e-300))
      }, numeric(1))
    )
  }))
  p_pc <- ggplot(pc_assoc, aes(x = PC, y = BatchVar, fill = neglog10p)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    scale_fill_gradient(low = "white", high = "#8338EC", name = "-log10(p)") +
    labs(title = "PC association with batch (Kruskal-Wallis)", x = "Principal component", y = NULL) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(out_dir, "batch_pc_association.png"), p_pc,
    width = 11, height = 5, dpi = 200, bg = "white"
  )

  message("\n=== Batch effect summary (Cramer's V, cluster vs batch) ===")
  for (i in seq_len(nrow(assoc_df))) {
    message(sprintf(
      "  %-18s V = %.3f  (p = %.2e)",
      assoc_df$BatchVar[i], assoc_df$cramers_v[i], assoc_df$p_value[i]
    ))
  }
  invisible(assoc_df)
}
