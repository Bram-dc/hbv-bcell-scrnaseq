# GSEA from a DEG table: rank genes, run clusterProfiler::GSEA, write table + plots.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(clusterProfiler)
  library(enrichplot)
  library(stringr)
})

run_gsea <- function(de, term2gene, out_dir, spec, label = spec$label) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # Drop ribosomal (RPL/RPS) and mitochondrial (MT-) genes so neither drives artifact enrichment (lab convention).
  rb_genes <- de$gene[grepl("^RPL|^RPS|^MT-", de$gene)]
  message(sprintf("  ribosomal/mitochondrial genes excluded from GSEA: %d", length(rb_genes)))

  gene_rank <- de |>
    filter(!gene %in% rb_genes) |>
    mutate(rank = sign(avg_log2FC) * -log10(p_val + 1e-300)) |>
    arrange(desc(rank))
  ranked <- setNames(gene_rank$rank, gene_rank$gene)

  res <- tryCatch(
    GSEA(
      geneList = ranked, TERM2GENE = term2gene,
      minGSSize = GSEA_MIN_GS, maxGSSize = GSEA_MAX_GS,
      pvalueCutoff = GSEA_PCUT, verbose = FALSE, seed = TRUE
    ),
    error = function(e) {
      message(sprintf("  GSEA failed: %s", e$message))
      NULL
    }
  )
  if (is.null(res)) {
    return(invisible(NULL))
  }

  gsea_df <- as.data.frame(res)
  message(sprintf("  GSEA combined: %d significant sets", nrow(gsea_df)))
  write.csv(gsea_df, file.path(out_dir, "gsea_combined.csv"), row.names = FALSE)
  if (nrow(gsea_df) == 0) {
    return(invisible(res))
  }

  # Truncate and de-duplicate long pathway descriptions for plotting.
  res_plot <- res
  trunc <- str_trunc(res_plot@result$Description, width = 60, ellipsis = "...")
  dupes <- ave(trunc, trunc, FUN = seq_along)
  res_plot@result$Description <- ifelse(dupes == "1", trunc, paste0(trunc, " [", dupes, "]"))

  n_terms <- min(25, nrow(gsea_df))
  plot_h <- max(8, n_terms * 0.35)
  plot_df <- head(res_plot@result |> arrange(p.adjust), n_terms)

  p_dot <- dotplot(res_plot, showCategory = n_terms, font.size = 8) +
    labs(
      title = sprintf("%s - GSEA combined", label),
      subtitle = "ImmuneSigDB + Reactome + GO:BP"
    ) +
    theme(
      axis.text.y = element_text(size = 8),
      plot.title = element_text(size = 11),
      plot.subtitle = element_text(size = 8, colour = "grey40")
    )
  ggsave(file.path(out_dir, "gsea_combined_dotplot.png"), p_dot,
    width = 13, height = plot_h, dpi = 200, bg = "white"
  )

  p_nes <- ggplot(
    plot_df |> arrange(NES),
    aes(x = NES, y = reorder(Description, NES), fill = NES > 0)
  ) +
    geom_col(width = 0.7) +
    geom_vline(xintercept = 0, linewidth = 0.4, colour = "grey30") +
    scale_fill_manual(
      values = c("TRUE" = "#E63946", "FALSE" = "#457B9D"),
      labels = c(
        "TRUE" = paste("Higher in", spec$ident1),
        "FALSE" = paste("Higher in", spec$ident2)
      ), name = NULL
    ) +
    labs(
      title = sprintf("%s - NES", label),
      subtitle = sprintf(
        "Positive = upregulated in %s  |  Negative = upregulated in %s",
        spec$ident1, spec$ident2
      ),
      x = "Normalised Enrichment Score", y = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.y = element_text(size = 8), legend.position = "bottom",
      plot.subtitle = element_text(size = 8, colour = "grey40")
    )
  ggsave(file.path(out_dir, "gsea_combined_nes.png"), p_nes,
    width = 11, height = plot_h, dpi = 200, bg = "white"
  )

  invisible(res)
}

# Drive GSEA off a DEG output tree, one gsea unit per de_all.csv.
run_gsea_from_deg <- function(spec, term2gene,
                              deg_root = file.path(OUTPUT_ROOT, spec$name, "deg"),
                              gsea_root = file.path(OUTPUT_ROOT, spec$name, "gsea")) {
  message(sprintf("\n== GSEA: %s ==", spec$label))
  units <- list.dirs(deg_root, recursive = FALSE, full.names = FALSE)
  for (u in units) {
    f <- file.path(deg_root, u, "de_all.csv")
    if (!file.exists(f)) next
    de <- read.csv(f, stringsAsFactors = FALSE)
    if (nrow(de) == 0 || !"gene" %in% names(de)) next
    message(sprintf("-- %s --", u))
    disp <- if (u == "whole") "whole" else cluster_name(sub("^cluster_", "", u), spec$antigen)
    run_gsea(de, term2gene, file.path(gsea_root, u), spec,
      label = paste0(spec$label, " [", disp, "]")
    )
  }
  invisible(gsea_root)
}
