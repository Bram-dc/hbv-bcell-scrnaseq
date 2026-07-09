# Per-cell module scoring (AddModuleScore) for curated MODULE_SETS programs, tested across conditions, whole and per cluster.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# Condition effect on a score vector within a cell subset; returns n1/n2, medians, delta, p.
# Batch as lm covariate only when >=2 shared batches, else Wilcoxon (batch-confounding guard).
.module_test <- function(score, cv, bv, ident1, ident2) {
  cv <- droplevels(cv)
  bv <- droplevels(bv)
  n1 <- sum(cv == ident1)
  n2 <- sum(cv == ident2)
  if (n1 < 3 || n2 < 3) {
    return(data.frame(n1 = n1, n2 = n2, median1 = NA, median2 = NA, delta = NA, p = NA))
  }
  tab <- table(cv, bv)
  shared <- if (nrow(tab) == 2) sum(colSums(tab > 0) == 2) else 0
  p <- tryCatch(
    if (shared >= 2) {
      coef(summary(lm(score ~ cv + bv)))[2, "Pr(>|t|)"]
    } else {
      wilcox.test(score[cv == ident1], score[cv == ident2])$p.value
    },
    error = function(e) NA_real_
  )
  m1 <- median(score[cv == ident1], na.rm = TRUE)
  m2 <- median(score[cv == ident2], na.rm = TRUE)
  data.frame(n1 = n1, n2 = n2, median1 = m1, median2 = m2, delta = m1 - m2, p = p)
}

run_modules <- function(obj, spec, module_sets = MODULE_SETS,
                        out_dir = file.path(OUTPUT_ROOT, spec$name, "modules")) {
  message(sprintf("\n== Module scores: %s ==", spec$label))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # Score each module once on the whole comparison object (genes present only).
  score_cols <- character(0)
  for (nm in names(module_sets)) {
    genes <- intersect(module_sets[[nm]], rownames(obj))
    if (length(genes) < 3) {
      message(sprintf("  module %-18s <3 of %d genes present - skipped", nm, length(module_sets[[nm]])))
      next
    }
    obj <- AddModuleScore(obj,
      features = list(genes), name = paste0("mod_", nm, "_"),
      seed = 42, search = FALSE
    )
    score_cols[nm] <- paste0("mod_", nm, "_1")
    message(sprintf("  module %-18s %d/%d genes", nm, length(genes), length(module_sets[[nm]])))
  }
  if (length(score_cols) == 0) {
    message("  no modules scored")
    return(invisible(NULL))
  }

  meta <- obj@meta.data
  cond <- factor(as.character(obj$condition), levels = c(spec$ident1, spec$ident2))
  batch <- factor(trimws(as.character(obj@meta.data[[BATCH_VAR]])))
  clust <- as.character(obj$seurat_clusters)

  # Whole + (optional) per-cluster index sets.
  units <- c(
    list(whole = seq_len(ncol(obj))),
    if (isTRUE(spec$per_cluster)) split(seq_len(ncol(obj)), clust) else NULL
  )

  rows <- list()
  for (u in names(units)) {
    idx <- units[[u]]
    for (nm in names(score_cols)) {
      res <- .module_test(
        meta[[score_cols[nm]]][idx], cond[idx], batch[idx],
        spec$ident1, spec$ident2
      )
      rows[[length(rows) + 1]] <- cbind(unit = u, module = nm, res)
    }
  }
  summ <- do.call(rbind, rows)
  summ$p_adj <- p.adjust(summ$p, method = "BH")
  write.csv(summ, file.path(out_dir, "module_scores_summary.csv"), row.names = FALSE)

  # Whole-comparison violins (one panel per module).
  plot_df <- data.frame(
    Condition = cond,
    setNames(as.data.frame(lapply(score_cols, function(c) meta[[c]])), names(score_cols)),
    check.names = FALSE
  )
  pl <- plot_df |> pivot_longer(-Condition, names_to = "module", values_to = "score")
  pal <- setNames(c("#E63946", "#457B9D"), c(spec$ident1, spec$ident2))
  p <- ggplot(pl, aes(Condition, score, fill = Condition)) +
    geom_violin(trim = TRUE, scale = "width", alpha = 0.7) +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, colour = "grey30") +
    facet_wrap(~module, scales = "free_y") +
    scale_fill_manual(values = pal) +
    labs(
      title = paste0(spec$label, " - targeted module scores"),
      subtitle = "Per-cell AddModuleScore; box = IQR/median", x = NULL, y = "Module score"
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none", strip.text = element_text(face = "bold"))
  ggsave(file.path(out_dir, "module_scores_violin.png"), p,
    width = 11, height = 8,
    dpi = 200, bg = "white"
  )

  # Per-cluster delta heatmap (median ident1 - ident2), * = BH-adj p < 0.05.
  if (isTRUE(spec$per_cluster)) {
    hm <- summ[summ$unit != "whole" & !is.na(summ$delta), ]
    if (nrow(hm) > 0) {
      hm$sig <- ifelse(!is.na(hm$p_adj) & hm$p_adj < 0.05, "*", "")
      lev_ids <- as.character(sort(unique(as.integer(hm$unit))))
      hm$unit <- factor(cluster_name(hm$unit, spec$antigen),
        levels = cluster_name(lev_ids, spec$antigen)
      )
      ph <- ggplot(hm, aes(unit, module, fill = delta)) +
        geom_tile(colour = "white", linewidth = 0.4) +
        geom_text(aes(label = sig), size = 5, vjust = 0.75) +
        scale_fill_gradient2(low = "#457B9D", mid = "white", high = "#E63946", midpoint = 0) +
        labs(
          title = sprintf("%s - module Δ (%s − %s) per subset", spec$label, spec$ident1, spec$ident2),
          subtitle = "* BH-adjusted p < 0.05", x = NULL, y = NULL, fill = "Δ median"
        ) +
        theme_minimal(base_size = 11) +
        theme(
          panel.grid = element_blank(),
          axis.text.x = element_text(angle = 20, hjust = 1)
        )
      ggsave(file.path(out_dir, "module_scores_per_cluster.png"), ph,
        width = 9, height = 5,
        dpi = 200, bg = "white"
      )
    }
  }
  message(sprintf("  scored %d modules -> %s", length(score_cols), out_dir))
  invisible(summ)
}
