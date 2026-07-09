# Reusable plotters: DEG volcano/barplot and cluster-EDA panels (UMAPs, composition, QC violins, marker dotplot, feature maps).

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(ggExtra)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(EnhancedVolcano)
})

.cluster_factor <- function(obj) {
  if (!is.null(obj@meta.data[["cluster_label"]])) {
    cl <- obj@meta.data[["cluster_label"]]
    return(factor(as.character(cl), levels = levels(cl)))
  }
  cl <- as.character(obj$seurat_clusters)
  factor(cl, levels = as.character(sort(unique(as.integer(cl)))))
}

# Ranked by the GSEA metric so volcano and barplot highlight the same genes.
top_deg_genes <- function(de, n = TOP_N) {
  scored <- de |>
    filter(p_val < DEG_P, avg_log2FC != 0) |>
    mutate(.score = sign(avg_log2FC) * -log10(p_val + 1e-300))
  up <- scored |>
    filter(avg_log2FC > 0) |>
    arrange(desc(.score), desc(avg_log2FC)) |>
    head(n)
  down <- scored |>
    filter(avg_log2FC < 0) |>
    arrange(.score, avg_log2FC) |>
    head(n)
  bind_rows(up, down) |> select(-.score)
}

volcano_plot <- function(de, spec, n1, n2, out_file) {
  top_labels <- top_deg_genes(de) |> pull(gene)

  p <- EnhancedVolcano(de,
    lab = de$gene, x = "avg_log2FC", y = "p_val",
    pCutoff = DEG_P, FCcutoff = DEG_LOGFC, selectLab = top_labels,
    pointSize = 2.0, labSize = 3.5, drawConnectors = TRUE,
    widthConnectors = 0.4, colConnectors = "grey50",
    col = c("grey70", "#A8DADC", "#457B9D", "#E63946"),
    title = spec$label,
    subtitle = sprintf(
      "%s: %d cells  |  %s: %d cells  |  p_val < %.2f  |  |log2FC| > %.2f",
      spec$ident1, n1, spec$ident2, n2, DEG_P, DEG_LOGFC
    ),
    xlab = sprintf("avg log2FC  (positive = higher in %s)", spec$ident1),
    legendLabels = c("NS", "|log2FC|", "p_val", "p_val & |log2FC|")
  )
  ggsave(out_file, p, width = 10, height = 8, dpi = 200, bg = "white")
}

de_barplot <- function(de_sig, spec, n1, n2, out_file) {
  top_genes <- top_deg_genes(de_sig)
  if (nrow(top_genes) == 0) {
    return(invisible(NULL))
  }

  p <- ggplot(top_genes, aes(
    x = avg_log2FC, y = reorder(gene, avg_log2FC),
    fill = avg_log2FC > 0
  )) +
    geom_bar(stat = "identity") +
    geom_vline(xintercept = 0, linewidth = 0.4, colour = "grey40") +
    scale_fill_manual(
      values = c("TRUE" = "#e63946", "FALSE" = "#457b9d"),
      labels = c(
        "TRUE" = paste("Higher in", spec$ident1),
        "FALSE" = paste("Higher in", spec$ident2)
      ),
      name = NULL
    ) +
    labs(
      title = spec$label,
      subtitle = sprintf(
        "%s: %d cells  |  %s: %d cells  |  p_val < %.2f",
        spec$ident1, n1, spec$ident2, n2, DEG_P
      ),
      x = sprintf("avg log2FC  (positive = higher in %s)", spec$ident1), y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")
  ggsave(out_file, p,
    width = 8, height = max(5, nrow(top_genes) * 0.28),
    dpi = 200, bg = "white"
  )
}

# Group-placement UMAP: cells coloured by their two conditions on the shared embedding.
umap_condition_plot <- function(obj, spec, out_file) {
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  df <- as.data.frame(Embeddings(obj, "umap"))
  colnames(df) <- c("UMAP_1", "UMAP_2")
  df$Condition <- factor(as.character(obj$condition), levels = c(spec$ident1, spec$ident2))
  n1 <- sum(obj$condition == spec$ident1)
  n2 <- sum(obj$condition == spec$ident2)
  pal <- setNames(c("#E63946", "#457B9D"), c(spec$ident1, spec$ident2))
  base <- ggplot(df, aes(UMAP_1, UMAP_2, colour = Condition)) +
    geom_point(alpha = 0.6, size = 1.2) +
    scale_colour_manual(values = pal) +
    labs(
      title = paste0(spec$label, " - group placement"),
      subtitle = sprintf(
        "%s (%d cells)  vs  %s (%d cells)  on the shared embedding",
        spec$ident1, n1, spec$ident2, n2
      ),
      colour = NULL
    ) +
    theme_minimal(base_size = 12) +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))
  # Per-group marginal density curves on both axes (matches the other UMAPs).
  p <- ggMarginal(base, type = "density", groupColour = TRUE, groupFill = TRUE)
  ggsave(out_file, p, width = 8, height = 6, dpi = 200, bg = "white")
}

# UMAP coloured by cluster (with marginal densities).
umap_cluster_plot <- function(obj, out_file, title = "Unsupervised clusters",
                              subtitle = NULL) {
  df <- as.data.frame(Embeddings(obj, "umap"))
  colnames(df) <- c("UMAP_1", "UMAP_2")
  df$Cluster <- .cluster_factor(obj)
  lev <- levels(df$Cluster)
  pal <- setNames(scales::hue_pal()(length(lev)), lev)
  if (is.null(subtitle)) {
    subtitle <- paste0(length(lev), " clusters,  ", nrow(df), " cells")
  }
  base <- ggplot(df, aes(UMAP_1, UMAP_2, colour = Cluster)) +
    geom_point(alpha = 0.6, size = 1.2) +
    scale_colour_manual(values = pal) +
    labs(title = title, subtitle = subtitle, colour = "Cluster") +
    theme_minimal(base_size = 12) +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))
  # ggMarginal wedges the legend between panel and density strip; re-attach it to the right of the whole marginal plot.
  legend <- cowplot::get_legend(base)
  marg <- ggMarginal(base + theme(legend.position = "none"),
    type = "density", groupColour = TRUE, groupFill = TRUE
  )
  p <- cowplot::plot_grid(marg, legend, nrow = 1, rel_widths = c(1, 0.32))
  ggsave(out_file, p, width = 9, height = 6, dpi = 200, bg = "white")
  invisible(pal)
}

# UMAP coloured by an arbitrary metadata column.
umap_meta_plot <- function(obj, col, out_file, title = NULL) {
  df <- as.data.frame(Embeddings(obj, "umap"))
  colnames(df) <- c("UMAP_1", "UMAP_2")
  df[[col]] <- trimws(as.character(obj@meta.data[[col]]))
  if (is.null(title)) title <- paste("UMAP by", col)
  p <- ggplot(df, aes(UMAP_1, UMAP_2, colour = .data[[col]])) +
    geom_point(alpha = 0.6, size = 1.2) +
    labs(
      title = title,
      subtitle = paste0(dplyr::n_distinct(df[[col]]), " levels,  ", nrow(df), " cells"),
      colour = col
    ) +
    theme_minimal(base_size = 12) +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))
  ggsave(out_file, p, width = 8, height = 6, dpi = 200, bg = "white")
}

# Stacked cluster composition per group.
composition_plot <- function(obj, out_file, title = "Cluster composition per group") {
  cf <- .cluster_factor(obj)
  lev <- levels(cf)
  pal <- setNames(scales::hue_pal()(length(lev)), lev)
  stats <- obj@meta.data |>
    group_by(Group = trimws(Group)) |>
    summarise(n_patients = n_distinct(Sample), n_cells = n(), .groups = "drop")
  df <- data.frame(
    Group = trimws(obj$Group),
    Cluster = as.character(cf)
  ) |>
    count(Group, Cluster) |>
    group_by(Group) |>
    mutate(Proportion = n / sum(n), Cluster = factor(Cluster, levels = lev)) |>
    ungroup() |>
    left_join(stats, by = "Group") |>
    mutate(Group = sprintf("%s\n%d patients, %d cells", Group, n_patients, n_cells))
  p <- ggplot(df, aes(x = Group, y = Proportion, fill = Cluster)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = pal) +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(title = title, x = NULL, y = "Proportion of cells", fill = "Cluster") +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major.x = element_blank()
    )
  ggsave(out_file, p,
    width = max(8, dplyr::n_distinct(df$Group) * 0.7), height = 6,
    dpi = 200, bg = "white"
  )
}

# Stacked sub-cluster composition collapsed to the three viral-activity states (IA/IC/FC).
composition_state_plot <- function(obj, out_file,
                                   title = "Sub-cluster composition per viral activity state") {
  cf <- .cluster_factor(obj)
  lev <- levels(cf)
  pal <- setNames(scales::hue_pal()(length(lev)), lev)
  gc <- trimws(obj$Group)
  state <- dplyr::case_when(
    gc == "IA" ~ "IA",
    gc %in% c("IC(high)", "IC(low)") ~ "IC",
    gc %in% FC_GROUPS ~ "FC",
    TRUE ~ NA_character_
  )
  base <- data.frame(state = state, Sample = obj$Sample, Cluster = as.character(cf)) |>
    dplyr::filter(!is.na(state))
  stats <- base |>
    dplyr::group_by(state) |>
    dplyr::summarise(n_patients = dplyr::n_distinct(Sample), n_cells = dplyr::n(), .groups = "drop")
  df <- base |>
    dplyr::count(state, Cluster) |>
    dplyr::group_by(state) |>
    dplyr::mutate(Proportion = n / sum(n), Cluster = factor(Cluster, levels = lev)) |>
    dplyr::ungroup() |>
    dplyr::left_join(stats, by = "state") |>
    dplyr::mutate(
      state = factor(state, levels = c("IA", "IC", "FC")),
      label = sprintf(
        "%s\n%d patients, %s cells", state, n_patients,
        formatC(n_cells, big.mark = ",", format = "d")
      )
    )
  lab_lev <- df |>
    dplyr::distinct(state, label) |>
    dplyr::arrange(state) |>
    dplyr::pull(label)
  df$label <- factor(df$label, levels = lab_lev)
  p <- ggplot(df, aes(x = label, y = Proportion, fill = Cluster)) +
    geom_bar(stat = "identity", width = 0.8) +
    scale_fill_manual(values = pal) +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(title = title, x = NULL, y = "Proportion of cells", fill = "Sub-cluster") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.major.x = element_blank())
  ggsave(out_file, p, width = 8, height = 5, dpi = 200, bg = "white")
}

# Cohort overview collapsed to the three viral-activity states: (a) patients per state, (b) cells per patient by sub-cluster.
cohort_state_plot <- function(obj, out_file) {
  gc <- trimws(obj$Group)
  state <- dplyr::case_when(
    gc == "IA" ~ "IA",
    gc %in% c("IC(high)", "IC(low)") ~ "IC",
    gc %in% FC_GROUPS ~ "FC",
    TRUE ~ NA_character_
  )
  cf <- .cluster_factor(obj)
  md <- data.frame(state = state, Sample = obj$Sample, Cluster = cf) |>
    dplyr::filter(!is.na(state)) |>
    dplyr::mutate(state = factor(state, levels = c("IA", "IC", "FC")))
  state_pal <- c("IA" = "#E63946", "IC" = "#F4A261", "FC" = "#457B9D")
  clev <- levels(cf)
  clus_pal <- setNames(scales::hue_pal()(length(clev)), clev)

  pat <- md |>
    dplyr::distinct(state, Sample) |>
    dplyr::count(state, name = "n")
  pa <- ggplot(pat, aes(state, n, fill = state)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = n), vjust = -0.3, size = 4) +
    scale_fill_manual(values = state_pal, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(title = "(a) Patients per state", x = NULL, y = "Patients") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.major.x = element_blank())

  cellc <- md |> dplyr::count(state, Sample, Cluster, name = "n")
  ord <- cellc |>
    dplyr::group_by(state, Sample) |>
    dplyr::summarise(total = sum(n), .groups = "drop") |>
    dplyr::arrange(state, dplyr::desc(total))
  cellc$Sample <- factor(cellc$Sample, levels = ord$Sample)
  pb <- ggplot(cellc, aes(Sample, n, fill = Cluster)) +
    geom_col(width = 0.9) +
    facet_grid(~state, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = clus_pal, name = "Sub-cluster") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(title = "(b) Cells per patient by sub-cluster", x = "Patient (one bar each)", y = "Cells") +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      panel.grid.major.x = element_blank()
    )

  p <- patchwork::wrap_plots(pa, pb, widths = c(1, 3))
  ggsave(out_file, p, width = 12, height = 4.6, dpi = 200, bg = "white")
}

# QC metrics per cluster (mt / rb / mRNA).
qc_violin_plot <- function(obj, out_file, title = "QC metrics per cluster") {
  cf <- .cluster_factor(obj)
  lev <- levels(cf)
  pal <- setNames(scales::hue_pal()(length(lev)), lev)
  df <- obj@meta.data |>
    mutate(
      Cluster = cf,
      percent.mRNA = nCount_RNA / (nCount_RNA + nCount_spike_ins) * 100
    ) |>
    select(Cluster, percent.mt, percent.rb, percent.mRNA) |>
    pivot_longer(c(percent.mt, percent.rb, percent.mRNA),
      names_to = "Metric", values_to = "Value"
    ) |>
    mutate(Metric = factor(Metric,
      levels = c("percent.mt", "percent.rb", "percent.mRNA"),
      labels = c("% Mitochondrial", "% Ribosomal", "% mRNA")
    ))
  p <- ggplot(df, aes(x = Cluster, y = Value, fill = Cluster)) +
    geom_violin(trim = TRUE, scale = "width", alpha = 0.75) +
    geom_boxplot(width = 0.08, fill = "white", outlier.shape = NA, colour = "grey30") +
    scale_fill_manual(values = pal) +
    facet_wrap(~Metric, scales = "free_y", ncol = 1) +
    labs(
      title = title, subtitle = "Violin = distribution, box = IQR/median",
      x = "Cluster", y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none", strip.text = element_text(face = "bold", size = 10))
  ggsave(out_file, p, width = 12, height = 10, dpi = 200, bg = "white")
}

# Marker dotplot: scaled mean expression + % expressing, per cluster.
marker_dotplot <- function(obj, markers, out_file,
                           title = "Marker gene expression per cluster",
                           subtitle = "Dot size = % cells expressing,  Colour = scaled mean expression") {
  markers <- markers[markers %in% rownames(obj)]
  expr <- as.matrix(GetAssayData(obj, layer = "data"))
  cf <- .cluster_factor(obj)
  lev <- levels(cf)
  df <- do.call(rbind, lapply(lev, function(cl) {
    cells <- colnames(obj)[as.character(cf) == cl]
    do.call(rbind, lapply(markers, function(g) {
      data.frame(
        Cluster = cl, Gene = g,
        mean_expr = mean(expr[g, cells]),
        pct_express = mean(expr[g, cells] > 0) * 100
      )
    }))
  })) |>
    group_by(Gene) |>
    mutate(scaled_expr = (mean_expr - min(mean_expr)) /
      (max(mean_expr) - min(mean_expr) + 1e-9)) |>
    ungroup() |>
    mutate(
      Cluster = factor(Cluster, levels = lev),
      Gene = factor(Gene, levels = rev(markers))
    )
  p <- ggplot(df, aes(x = Cluster, y = Gene)) +
    geom_point(aes(size = pct_express, colour = scaled_expr)) +
    scale_colour_gradient(low = "grey92", high = "#b23a48", name = "Scaled\nmean expr") +
    scale_size_continuous(range = c(0.5, 8), name = "% expressing") +
    labs(title = title, subtitle = subtitle, x = "Cluster", y = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.y = element_text(size = 10),
      panel.grid.major = element_line(colour = "grey93"), legend.position = "right"
    )
  ggsave(out_file, p, width = 12, height = 8, dpi = 200, bg = "white")
}

# Violin of key markers per cluster.
violin_markers_plot <- function(obj, markers, out_file,
                                title = "B cell marker expression per cluster") {
  markers <- markers[markers %in% rownames(obj)]
  expr <- as.matrix(GetAssayData(obj, layer = "data"))
  cf <- .cluster_factor(obj)
  lev <- levels(cf)
  pal <- setNames(scales::hue_pal()(length(lev)), lev)
  df <- as.data.frame(t(expr[markers, , drop = FALSE])) |>
    tibble::rownames_to_column("cell") |>
    mutate(Cluster = cf) |>
    pivot_longer(all_of(markers), names_to = "Gene", values_to = "Expression") |>
    mutate(Gene = factor(Gene, levels = markers))
  p <- ggplot(df, aes(x = Cluster, y = Expression, fill = Cluster)) +
    geom_violin(trim = TRUE, scale = "width", alpha = 0.75) +
    geom_boxplot(width = 0.08, fill = "white", outlier.shape = NA, colour = "grey30") +
    scale_fill_manual(values = pal) +
    facet_wrap(~Gene, scales = "free_y", ncol = 2) +
    labs(
      title = title, subtitle = "Violin = distribution, box = IQR/median",
      x = "Cluster", y = "Normalised expression"
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none", strip.text = element_text(face = "bold", size = 10))
  ggsave(out_file, p,
    width = 14, height = ceiling(length(markers) / 2) * 3,
    dpi = 200, bg = "white"
  )
}
