# Per-cell ssGSEA (GSVA) scoring and condition test for a comparison's two conditions; lm ~ condition + batch, Wilcoxon fallback when batch is confounded.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(GSVA)
})

.ssgsea_boxplot <- function(score_long, summary_tbl, title, out_file) {
  top <- summary_tbl |>
    filter(!is.na(p_adj)) |>
    slice_head(n = 8) |>
    pull(pathway)
  if (length(top) == 0) {
    top <- summary_tbl |>
      slice_max(order_by = abs(delta_median), n = 8, with_ties = FALSE) |>
      pull(pathway)
  }
  df <- score_long |>
    filter(pathway %in% top) |>
    mutate(pathway = factor(pathway, levels = rev(top)))
  p <- ggplot(df, aes(x = Condition, y = score, fill = Condition)) +
    geom_violin(trim = TRUE, scale = "width", alpha = 0.7) +
    geom_boxplot(width = 0.08, fill = "white", outlier.shape = NA, colour = "grey30") +
    facet_wrap(~pathway, scales = "free_y", ncol = 2) +
    scale_fill_brewer(palette = "Set2") +
    labs(
      title = paste0("ssGSEA scores: ", title),
      subtitle = "Top pathways ranked by BH-adjusted condition p-value",
      x = NULL, y = "ssGSEA score"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none", strip.text = element_text(face = "bold"))
  ggsave(out_file, p, width = 12, height = 9)
}

.ssgsea_heatmap <- function(score_long, title, out_file) {
  top <- score_long |>
    group_by_at("pathway") |>
    summarise(variance = var(score, na.rm = TRUE), .groups = "drop") |>
    slice_max(order_by = variance, n = 12, with_ties = FALSE) |>
    pull(pathway)
  df <- score_long |>
    filter(pathway %in% top) |>
    mutate(pathway = factor(pathway, levels = rev(top)))
  p <- ggplot(df, aes(x = cell, y = pathway, fill = score)) +
    geom_tile() +
    facet_grid(cols = vars(Condition), scales = "free_x", space = "free_x") +
    scale_fill_gradient2(low = "#30638e", mid = "#f4f4f4", high = "#b23a48", midpoint = 0) +
    labs(
      title = paste0("ssGSEA heatmap: ", title),
      subtitle = "Top variable pathways; cells grouped by condition",
      x = NULL, y = NULL, fill = "Score"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold"),
      panel.spacing.x = grid::unit(6, "pt")
    )
  ggsave(out_file, p, width = 14, height = 8)
}

# ssGSEA for one collection on a two-condition comparison object.
.ssgsea_collection <- function(obj, expr, coll, spec, out_dir) {
  coll_dir <- file.path(out_dir, coll$slug)
  dir.create(coll_dir, recursive = TRUE, showWarnings = FALSE)

  genes_in_sets <- unique(unlist(coll$gene_sets, use.names = FALSE))
  common <- intersect(rownames(expr), genes_in_sets)
  if (length(common) < 10) {
    message("Skipping ", coll$label, " - too few overlapping genes")
    return(invisible(NULL))
  }
  sets <- lapply(coll$gene_sets, function(s) intersect(s, common))
  sets <- sets[lengths(sets) >= 10]
  expr_sub <- expr[common, , drop = FALSE]

  scores <- gsva(ssgseaParam(
    exprData = expr_sub, geneSets = sets,
    minSize = 10, maxSize = Inf, normalize = TRUE
  ), verbose = FALSE)

  slug <- paste0(spec$name, "_", coll$slug)
  saveRDS(scores, file.path(coll_dir, paste0("ssgsea_scores_", slug, ".Rds")))
  write.csv(as.data.frame(scores), file.path(coll_dir, paste0("ssgsea_scores_", slug, ".csv")))

  cell_meta <- tibble(
    cell = colnames(obj),
    Condition = factor(as.character(obj$condition), levels = c(spec$ident1, spec$ident2)),
    Batch = trimws(as.character(obj@meta.data[[BATCH_VAR]]))
  )

  score_long <- as.data.frame(scores) |>
    tibble::rownames_to_column("pathway") |>
    pivot_longer(-all_of("pathway"), names_to = "cell", values_to = "score") |>
    left_join(cell_meta, by = "cell")

  summary_tbl <- score_long |>
    group_by_at("pathway") |>
    group_modify(~ {
      cv <- .x$Condition
      sc <- .x$score
      bv <- droplevels(factor(.x$Batch))
      n1 <- sum(cv == spec$ident1)
      n2 <- sum(cv == spec$ident2)
      tab <- table(cv, bv)
      shared <- if (nrow(tab) == 2) sum(colSums(tab > 0) == 2) else 0
      p <- if (n1 >= 2 && n2 >= 2) {
        tryCatch(
          if (shared >= 2) {
            coef(summary(lm(sc ~ cv + bv)))[2, "Pr(>|t|)"]
          } else {
            wilcox.test(sc[cv == spec$ident1], sc[cv == spec$ident2])$p.value
          },
          error = function(e) NA_real_
        )
      } else {
        NA_real_
      }
      tibble(
        n1 = n1, n2 = n2,
        median_1 = median(sc[cv == spec$ident1], na.rm = TRUE),
        median_2 = median(sc[cv == spec$ident2], na.rm = TRUE),
        p_value = p
      )
    }) |>
    ungroup() |>
    mutate(
      delta_median = median_1 - median_2,
      p_adj = p.adjust(p_value, method = "BH")
    ) |>
    arrange(p_adj, desc(abs(delta_median)))

  write.csv(summary_tbl, file.path(coll_dir, paste0("ssgsea_summary_", slug, ".csv")),
    row.names = FALSE
  )

  title <- paste0(spec$label, " (", coll$label, ")")
  .ssgsea_boxplot(
    score_long, summary_tbl, title,
    file.path(coll_dir, paste0("ssgsea_boxplot_", slug, ".png"))
  )
  .ssgsea_heatmap(
    score_long, title,
    file.path(coll_dir, paste0("ssgsea_heatmap_", slug, ".png"))
  )
  message("ssGSEA done: ", spec$label, " [", coll$label, "]")
  invisible(score_long)
}

run_ssgsea <- function(obj, spec, collections,
                       out_dir = file.path(OUTPUT_ROOT, spec$name, "ssgsea")) {
  message(sprintf("\n== ssGSEA: %s ==", spec$label))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (ncol(obj) == 0) stop("No cells for ssGSEA: ", spec$label)

  obj <- NormalizeData(obj, verbose = FALSE)
  expr <- GetAssayData(obj, layer = "data")
  expr <- expr[rowSums(expr > 0) > 0, , drop = FALSE]

  for (coll in collections) {
    .ssgsea_collection(obj, expr, coll, spec, out_dir)
  }
  invisible(out_dir)
}
