# Differential expression: writes DE tables + volcano + barplot per comparison.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

# LR batch covariate only when >=2 batches have both groups; else confounded, so skipped.
find_markers_batch <- function(object, ident.1, ident.2, batch_var = BATCH_VAR, ...) {
  b <- droplevels(factor(trimws(as.character(object@meta.data[[batch_var]]))))
  ids <- as.character(Idents(object))
  keep <- ids %in% c(ident.1, ident.2)
  tab <- table(ids[keep], b[keep])
  shared <- if (nrow(tab) == 2) sum(colSums(tab > 0) == 2) else 0
  if (shared >= 2) {
    object@meta.data[[batch_var]] <- b
    message(sprintf(
      "    [batch] %s vs %s: adjusting for %s (%d shared batches)",
      ident.1, ident.2, batch_var, shared
    ))
    FindMarkers(object,
      ident.1 = ident.1, ident.2 = ident.2,
      test.use = "LR", latent.vars = batch_var, ...
    )
  } else {
    message(sprintf(
      "    [batch] %s vs %s: confounded with %s (%d shared batches) - covariate skipped",
      ident.1, ident.2, batch_var, shared
    ))
    FindMarkers(object, ident.1 = ident.1, ident.2 = ident.2, ...)
  }
}

# One DE table for a (sub)object on its $condition idents; writes tables + both plots.
deg_one <- function(obj, spec, out_dir) {
  Idents(obj) <- obj$condition
  n1 <- sum(obj$condition == spec$ident1)
  n2 <- sum(obj$condition == spec$ident2)
  message(sprintf("  %s: %d cells  |  %s: %d cells", spec$ident1, n1, spec$ident2, n2))

  de <- find_markers_batch(obj,
    ident.1 = spec$ident1, ident.2 = spec$ident2,
    logfc.threshold = DEG_LOGFC, min.pct = DEG_MIN_PCT, verbose = FALSE
  )
  de$gene <- rownames(de)
  de_sig <- de |>
    filter(p_val < DEG_P) |>
    arrange(desc(avg_log2FC))

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(de, file.path(out_dir, "de_all.csv"), row.names = FALSE)
  write.csv(de_sig, file.path(out_dir, "de_significant.csv"), row.names = FALSE)

  volcano_plot(de, spec, n1, n2, file.path(out_dir, "volcano.png"))
  if (nrow(de_sig) > 0) {
    de_barplot(de_sig, spec, n1, n2, file.path(out_dir, "de_barplot.png"))
  }
  message(sprintf("  significant DE genes: %d", nrow(de_sig)))
  invisible(de)
}

# Run DEG for a comparison: whole object, plus each cluster if per_cluster.
run_deg <- function(obj, spec, out_dir = file.path(OUTPUT_ROOT, spec$name, "deg")) {
  message(sprintf("\n== DEG: %s ==", spec$label))
  deg_one(obj, spec, file.path(out_dir, "whole"))
  if (isTRUE(spec$per_cluster)) {
    for (cl in cluster_ids(obj)) {
      message(sprintf("-- cluster %s --", cl))
      deg_one(subset_cluster(obj, cl), spec, file.path(out_dir, paste0("cluster_", cl)))
    }
  }
  invisible(out_dir)
}
