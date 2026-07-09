# Cluster marker discovery (FindAllMarkers) and top-DE-gene dotplots.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(msigdbr)
})

# Internal: top-N DE-gene dotplot restricted to a gene universe.
.de_dotplot <- function(obj_de, cluster_de, expr, universe, lnc, coding, label, out_file) {
  lev <- as.character(sort(as.integer(unique(as.character(obj_de$seurat_clusters)))))
  lnc_de <- intersect(lnc, cluster_de$gene)
  genes <- if (is.null(universe)) {
    union(coding, lnc_de)
  } else {
    union(intersect(coding, universe), lnc_de)
  }
  de_filt <- cluster_de |>
    filter(gene %in% genes, pct.1 >= MIN_PCT_EXPRESS / 100, p_val_adj < 0.05)
  top <- de_filt |>
    group_by(cluster) |>
    slice_max(avg_log2FC, n = TOP_DE_N, with_ties = FALSE) |>
    ungroup()
  sel <- unique(top$gene)
  message(sprintf(
    "%s: %d top DE genes across %d clusters", label, length(sel),
    dplyr::n_distinct(top$cluster)
  ))
  if (length(sel) == 0) {
    return(invisible(NULL))
  }

  df <- do.call(rbind, lapply(lev, function(cl) {
    cells <- colnames(obj_de)[as.character(obj_de$seurat_clusters) == cl]
    do.call(rbind, lapply(sel, function(g) {
      data.frame(
        Cluster = cl, Gene = g,
        mean_expr = mean(expr[g, cells]), pct_express = mean(expr[g, cells] > 0) * 100
      )
    }))
  })) |>
    group_by(Gene) |>
    mutate(scaled_expr = (mean_expr - min(mean_expr)) /
      (max(mean_expr) - min(mean_expr) + 1e-9)) |>
    ungroup() |>
    mutate(Cluster = factor(Cluster, levels = lev), Gene = factor(Gene, levels = rev(sel)))
  p <- ggplot(df, aes(x = Cluster, y = Gene)) +
    geom_point(aes(size = pct_express, colour = scaled_expr)) +
    scale_colour_gradient(low = "grey92", high = "#b23a48", name = "Scaled\nmean expr") +
    scale_size_continuous(range = c(0.5, 8), name = "% expressing") +
    labs(
      title = paste0("Top ", TOP_DE_N, " DE genes per cluster (", label, ")"),
      subtitle = "Dot size = % expressing  |  Colour = scaled mean expression",
      x = "Seurat cluster", y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.y = element_text(size = 9),
      panel.grid.major = element_line(colour = "grey93"), legend.position = "right"
    )
  ggsave(out_file, p, width = 12, height = max(8, length(sel) * 0.32), dpi = 200, bg = "white")
  write.csv(de_filt, sub("\\.png$", ".csv", out_file), row.names = FALSE)
}

# FindAllMarkers (raw + filtered) plus the three DE dotplots.
cluster_markers <- function(obj, out_dir, exclude_clusters = character(0)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  expr <- as.matrix(GetAssayData(obj, layer = "data"))

  Idents(obj) <- obj$seurat_clusters
  de_raw <- FindAllMarkers(obj, only.pos = TRUE, logfc.threshold = DEG_LOGFC, verbose = FALSE)
  write.csv(de_raw, file.path(out_dir, "de_all_clusters_raw.csv"), row.names = FALSE)
  message(sprintf("Raw DE results exported: %d rows", nrow(de_raw)))

  obj_de <- if (length(exclude_clusters) > 0) {
    subset(obj, cells = colnames(obj)[!as.character(obj$seurat_clusters) %in% exclude_clusters])
  } else {
    obj
  }
  Idents(obj_de) <- obj_de$seurat_clusters
  cluster_de <- FindAllMarkers(obj_de, only.pos = TRUE, logfc.threshold = DEG_LOGFC, verbose = FALSE)
  if (nrow(cluster_de) == 0) {
    return(invisible(NULL))
  }

  coding <- protein_coding_genes(unique(cluster_de$gene))
  lnc <- lnc_genes()

  # GO immune universe
  go_terms <- c(
    "GO:0006955", "GO:0002376", "GO:0042113", "GO:0030183",
    "GO:0045321", "GO:0002250", "GO:0019724", "GO:0050871"
  )
  go_annot <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
    keys = go_terms, columns = "SYMBOL", keytype = "GOALL"
  ))
  go_immune <- unique(go_annot$SYMBOL[!is.na(go_annot$SYMBOL)])

  # MSigDB C7/C2 B cell universe
  c7 <- msigdbr(species = "Homo sapiens", collection = "C7", subcollection = "IMMUNESIGDB")
  c2 <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME")
  msigdb_bcell <- union(
    c7$gene_symbol[grepl("BCELL|B_CELL", c7$gs_name)],
    c2$gene_symbol[grepl("B_CELL|_BCR_", c2$gs_name)]
  )

  .de_dotplot(
    obj_de, cluster_de, expr, go_immune, lnc, coding, "GO immune",
    file.path(out_dir, "dotplot_top_de_go_immune.png")
  )
  .de_dotplot(
    obj_de, cluster_de, expr, msigdb_bcell, lnc, coding, "MSigDB C7/C2",
    file.path(out_dir, "dotplot_top_de_msigdb.png")
  )
  .de_dotplot(
    obj_de, cluster_de, expr, NULL, lnc, coding, "top significant",
    file.path(out_dir, "dotplot_top_de_significant.png")
  )
  invisible(cluster_de)
}
