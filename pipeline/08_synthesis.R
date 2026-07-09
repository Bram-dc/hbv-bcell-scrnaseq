# Cross-comparison synthesis of IA-vs-FC against IC-vs-FC, whole HBcAg population and per shared cluster. -> output/synthesis/ia_ic_vs_fc/

source("config/config.R")
source("config/comparisons.R")
invisible(lapply(list.files("R", "[.]R$", full.names = TRUE), source))

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ggrepel)
})

IA <- "hbcag_ia_vs_fc"
IC <- "hbcag_ic_vs_fc"
OUT_DIR <- file.path(OUTPUT_ROOT, "synthesis", "ia_ic_vs_fc")
CLUSTERS <- c("0", "1", "2")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

read_de <- function(cmp, subdir) {
  f <- file.path(OUTPUT_ROOT, cmp, "deg", subdir, "de_all.csv")
  if (!file.exists(f)) {
    return(NULL)
  }
  read.csv(f, stringsAsFactors = FALSE)
}
read_gsea <- function(cmp, subdir) {
  f <- file.path(OUTPUT_ROOT, cmp, "gsea", subdir, "gsea_combined.csv")
  if (!file.exists(f) || file.size(f) < 5) {
    return(NULL)
  }
  read.csv(f, stringsAsFactors = FALSE)
}

cat_cols <- c(
  "Shared (same direction)" = "#2A9D8F", "Divergent (opposite)" = "#E63946",
  "IA-specific" = "#457B9D", "IC-specific" = "#F4A261"
)

# Compare IA-vs-FC against IC-vs-FC at one level; write compare tables and scatter plots, return category counts and Pearson concordance r.
process_level <- function(subdir, suffix, title_label) {
  ia <- read_de(IA, subdir)
  ic <- read_de(IC, subdir)
  if (is.null(ia) || is.null(ic)) {
    message(sprintf("-- %s --  missing DE table, skipping", suffix))
    return(NULL)
  }
  message(sprintf("-- %s --", suffix))

  sig <- function(d) d$p_val < DEG_P & abs(d$avg_log2FC) > DEG_LOGFC
  ia$sig <- sig(ia)
  ic$sig <- sig(ic)
  merged <- full_join(
    ia |> select(gene, ia_log2FC = avg_log2FC, ia_p = p_val, ia_sig = sig),
    ic |> select(gene, ic_log2FC = avg_log2FC, ic_p = p_val, ic_sig = sig),
    by = "gene"
  ) |>
    mutate(
      ia_sig = ifelse(is.na(ia_sig), FALSE, ia_sig),
      ic_sig = ifelse(is.na(ic_sig), FALSE, ic_sig),
      category = case_when(
        ia_sig & ic_sig & sign(ia_log2FC) == sign(ic_log2FC) ~ "Shared (same direction)",
        ia_sig & ic_sig & sign(ia_log2FC) != sign(ic_log2FC) ~ "Divergent (opposite)",
        ia_sig & !ic_sig ~ "IA-specific",
        !ia_sig & ic_sig ~ "IC-specific",
        TRUE ~ "NS in both"
      )
    )
  write.csv(merged, file.path(OUT_DIR, sprintf("de_compare_%s.csv", suffix)), row.names = FALSE)
  counts <- merged |>
    filter(category != "NS in both") |>
    count(category)

  de_both <- merged |> filter((ia_sig | ic_sig) & !is.na(ia_log2FC) & !is.na(ic_log2FC))
  de_n <- nrow(de_both)
  de_r <- if (de_n > 2) cor(de_both$ia_log2FC, de_both$ic_log2FC) else NA_real_
  message(sprintf(
    "  DE categories: %s | gene r=%s (n=%d)",
    paste(sprintf("%s=%d", counts$category, counts$n), collapse = ", "),
    ifelse(is.na(de_r), "NA", sprintf("%.3f", de_r)), de_n
  ))

  plot_df <- merged |>
    filter(ia_sig | ic_sig) |>
    mutate(
      ia_log2FC = ifelse(is.na(ia_log2FC), 0, ia_log2FC),
      ic_log2FC = ifelse(is.na(ic_log2FC), 0, ic_log2FC)
    )
  lab_df <- plot_df |>
    filter(category %in% c("Shared (same direction)", "Divergent (opposite)")) |>
    mutate(mag = abs(ia_log2FC) + abs(ic_log2FC)) |>
    arrange(desc(mag)) |>
    group_by(category) |>
    slice_head(n = 12) |>
    ungroup()
  rng <- range(c(plot_df$ia_log2FC, plot_df$ic_log2FC), na.rm = TRUE)
  de_sub <- if (is.na(de_r)) {
    "Each point = a gene significant in at least one contrast"
  } else {
    sprintf("Pearson r = %.2f over %d genes significant in >=1 contrast", de_r, de_n)
  }
  p <- ggplot(plot_df, aes(ia_log2FC, ic_log2FC, colour = category)) +
    geom_hline(yintercept = 0, colour = "grey80") +
    geom_vline(xintercept = 0, colour = "grey80") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
    geom_point(alpha = 0.7, size = 1.8) +
    ggrepel::geom_text_repel(
      data = lab_df, aes(label = gene), size = 3, max.overlaps = 20,
      show.legend = FALSE
    ) +
    scale_colour_manual(values = cat_cols, name = NULL) +
    coord_equal(xlim = rng, ylim = rng) +
    labs(
      title = sprintf("%s - IA vs IC differential response (vs FC)", title_label),
      subtitle = de_sub,
      x = "avg log2FC  (IA vs FC)", y = "avg log2FC  (IC vs FC)"
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")
  ggsave(file.path(OUT_DIR, sprintf("de_scatter_%s.png", suffix)), p,
    width = 8, height = 8, dpi = 200, bg = "white"
  )

  gcounts <- NULL
  gsea_n <- 0L
  gsea_r <- NA_real_
  ia_g <- read_gsea(IA, subdir)
  ic_g <- read_gsea(IC, subdir)
  if (!is.null(ia_g) && !is.null(ic_g) && nrow(ia_g) > 0 && nrow(ic_g) > 0) {
    gmerged <- full_join(
      ia_g |> select(ID, Description, ia_NES = NES, ia_padj = p.adjust),
      ic_g |> select(ID, ic_NES = NES, ic_padj = p.adjust),
      by = "ID"
    ) |>
      mutate(gcategory = case_when(
        !is.na(ia_NES) & !is.na(ic_NES) & sign(ia_NES) == sign(ic_NES) ~ "Shared (same direction)",
        !is.na(ia_NES) & !is.na(ic_NES) & sign(ia_NES) != sign(ic_NES) ~ "Divergent (opposite)",
        !is.na(ia_NES) ~ "IA-specific", TRUE ~ "IC-specific"
      ))
    write.csv(gmerged, file.path(OUT_DIR, sprintf("gsea_compare_%s.csv", suffix)), row.names = FALSE)
    gcounts <- gmerged |> count(gcategory)

    shared_g <- gmerged |> filter(gcategory %in% c("Shared (same direction)", "Divergent (opposite)"))
    gsea_n <- nrow(shared_g)
    gsea_r <- if (gsea_n > 2) cor(shared_g$ia_NES, shared_g$ic_NES) else NA_real_
    message(sprintf(
      "  GSEA categories: %s | NES r=%s (n=%d)",
      paste(sprintf("%s=%d", gcounts$gcategory, gcounts$n), collapse = ", "),
      ifelse(is.na(gsea_r), "NA", sprintf("%.3f", gsea_r)), gsea_n
    ))
    if (gsea_n > 0) {
      lab_g <- shared_g |>
        mutate(mag = abs(ia_NES) + abs(ic_NES)) |>
        arrange(desc(mag)) |>
        slice_head(n = 15) |>
        mutate(short = substr(Description, 1, 45))
      grng <- range(c(shared_g$ia_NES, shared_g$ic_NES), na.rm = TRUE)
      g_sub <- if (is.na(gsea_r)) {
        NULL
      } else {
        sprintf("Pearson r = %.2f over %d shared pathways", gsea_r, gsea_n)
      }
      pg <- ggplot(shared_g, aes(ia_NES, ic_NES, colour = gcategory)) +
        geom_hline(yintercept = 0, colour = "grey80") +
        geom_vline(xintercept = 0, colour = "grey80") +
        geom_abline(slope = 1, linetype = "dashed", colour = "grey60") +
        geom_point(alpha = 0.8, size = 2) +
        ggrepel::geom_text_repel(
          data = lab_g, aes(label = short), size = 2.6,
          max.overlaps = 20, show.legend = FALSE
        ) +
        scale_colour_manual(values = c(
          "Shared (same direction)" = "#2A9D8F",
          "Divergent (opposite)" = "#E63946"
        ), name = NULL) +
        coord_equal(xlim = grng, ylim = grng) +
        labs(
          title = sprintf("%s - shared GSEA pathways (NES)", title_label),
          subtitle = g_sub,
          x = "NES (IA vs FC)", y = "NES (IC vs FC)"
        ) +
        theme_minimal(base_size = 11) +
        theme(legend.position = "bottom")
      ggsave(file.path(OUT_DIR, sprintf("gsea_scatter_%s.png", suffix)), pg,
        width = 9, height = 9, dpi = 200, bg = "white"
      )
    }
  } else {
    message("  GSEA: one side empty - pathway comparison skipped")
  }

  list(
    de_counts = counts,
    gsea_counts = gcounts,
    r = data.frame(
      level = suffix, de_r = de_r, de_n = de_n,
      gsea_r = gsea_r, gsea_n = gsea_n
    )
  )
}

results <- list()
results[["whole"]] <- process_level("whole", "whole", "Whole population")
for (cl in CLUSTERS) {
  results[[paste0("cluster_", cl)]] <- process_level(
    paste0("cluster_", cl), paste0("cluster_", cl), cluster_name(cl, "HBcAg")
  )
}

de_summary <- list()
gsea_summary <- list()
for (cl in CLUSTERS) {
  res <- results[[paste0("cluster_", cl)]]
  if (is.null(res)) next
  if (!is.null(res$de_counts) && nrow(res$de_counts) > 0) {
    de_summary[[cl]] <- res$de_counts |> mutate(cluster = cl)
  }
  if (!is.null(res$gsea_counts) && nrow(res$gsea_counts) > 0) {
    gsea_summary[[cl]] <- res$gsea_counts |> mutate(cluster = cl)
  }
}

if (length(de_summary) > 0) {
  de_tab <- bind_rows(de_summary) |> pivot_wider(names_from = category, values_from = n, values_fill = 0)
  write.csv(de_tab, file.path(OUT_DIR, "de_category_summary.csv"), row.names = FALSE)
  cat("\n=== DE gene categories per cluster (IA vs IC, both vs FC) ===\n")
  print(as.data.frame(de_tab))
}
if (length(gsea_summary) > 0) {
  gsea_tab <- bind_rows(gsea_summary) |> pivot_wider(names_from = gcategory, values_from = n, values_fill = 0)
  write.csv(gsea_tab, file.path(OUT_DIR, "gsea_category_summary.csv"), row.names = FALSE)
  cat("\n=== GSEA pathway categories per cluster ===\n")
  print(as.data.frame(gsea_tab))
}

r_tab <- bind_rows(lapply(results, function(x) if (!is.null(x)) x$r))
write.csv(r_tab, file.path(OUT_DIR, "concordance_r.csv"), row.names = FALSE)
cat("\n=== Concordance (Pearson r) IA-vs-FC vs IC-vs-FC ===\n")
print(r_tab)

message(sprintf("\nSynthesis outputs saved to: %s", OUT_DIR))
