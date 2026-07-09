# SCENIC regulon inference and differential regulon activity, gated behind RUN_SCENIC and driven per comparison/cluster by pipeline/07_scenic.R.

# Differential regulon activity (Wilcoxon on AUC) between the two conditions.
scenic_differential <- function(run_dir, spec, out_dir) {
  suppressPackageStartupMessages({
    library(AUCell)
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(pheatmap)
    library(RColorBrewer)
  })
  run_dir <- normalizePath(run_dir, mustWork = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_dir <- normalizePath(out_dir, mustWork = FALSE)

  auc_file <- file.path(run_dir, "int", "3.4_regulonAUC.Rds")
  if (!file.exists(auc_file)) {
    stop("No cached regulon AUC at ", auc_file, " - run SCENIC inference first.")
  }
  regulonAUC <- readRDS(auc_file)
  cellInfo <- readRDS(file.path(run_dir, "int", "cellInfo.Rds"))
  auc_mat <- getAUC(regulonAUC)
  shared <- intersect(colnames(auc_mat), rownames(cellInfo))
  auc_mat <- auc_mat[, shared]
  cellInfo <- cellInfo[shared, ]

  g1 <- rownames(cellInfo)[cellInfo$Condition == spec$ident1]
  g2 <- rownames(cellInfo)[cellInfo$Condition == spec$ident2]
  message(sprintf(
    "%s: %d cells | %s: %d cells | regulons: %d",
    spec$ident1, length(g1), spec$ident2, length(g2), nrow(auc_mat)
  ))
  if (length(g1) < 3 || length(g2) < 3) {
    stop(sprintf(
      "Condition group too small: %s=%d, %s=%d - cache/spec mismatch?",
      spec$ident1, length(g1), spec$ident2, length(g2)
    ))
  }

  res <- do.call(rbind, lapply(rownames(auc_mat), function(reg) {
    x1 <- auc_mat[reg, g1]
    x2 <- auc_mat[reg, g2]
    wt <- wilcox.test(x1, x2, exact = FALSE)
    data.frame(
      p_value = wt$p.value, W = unname(wt$statistic),
      effect_rbc = as.numeric(1 - (2 * wt$statistic) / (length(x1) * length(x2))),
      mean_1 = mean(x1), mean_2 = mean(x2),
      log2FC = log2((mean(x1) + 1e-6) / (mean(x2) + 1e-6))
    )
  }))
  rownames(res) <- rownames(auc_mat)
  res$regulon <- rownames(res)
  res$p_adj <- p.adjust(res$p_value, method = "BH")
  res$sig <- res$p_adj < 0.05
  res <- res[order(res$p_adj), ]
  write.csv(res, file.path(out_dir, "wilcoxon_regulons.csv"), row.names = FALSE)
  message(sprintf("Significant regulons (FDR<0.05): %d", sum(res$sig, na.rm = TRUE)))

  tryCatch(
    {
      res$label_short <- ifelse(res$p_adj < 0.05, gsub(" \\(.*\\)", "", res$regulon), NA)
      idx <- which(!is.na(res$label_short))
      if (length(idx) > 15) res$label_short[setdiff(idx, head(idx, 15))] <- NA
      volcano <- ggplot(res, aes(x = log2FC, y = -log10(p_adj))) +
        geom_point(aes(color = sig, size = sig), alpha = 0.7) +
        scale_color_manual(
          values = c("FALSE" = "grey70", "TRUE" = "#c0392b"),
          labels = c("ns", "FDR < 0.05")
        ) +
        scale_size_manual(values = c("FALSE" = 1.2, "TRUE" = 2.2), guide = "none") +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
        ggrepel::geom_text_repel(aes(label = label_short), size = 3, max.overlaps = 20) +
        labs(
          title = sprintf("SCENIC regulon activity: %s", spec$label),
          x = sprintf("log2 fold change (%s / %s)", spec$ident1, spec$ident2),
          y = expression(-log[10](FDR)), color = NULL
        ) +
        theme_bw(base_size = 12) +
        theme(legend.position = "top")
      ggsave(file.path(out_dir, "volcano_regulons.pdf"), volcano, width = 8, height = 6)
    },
    error = function(e) message("volcano plot failed (non-fatal): ", conditionMessage(e))
  )

  tryCatch(
    {
      sig_regulons <- res$regulon[which(res$sig)]
      if (length(sig_regulons) == 0) sig_regulons <- head(res$regulon, 30)
      # res is ordered by p_adj, so head() keeps the top hits; cap the count to keep the heatmap legible.
      n_show <- min(length(sig_regulons), 40)
      sig_regulons <- head(sig_regulons, n_show)
      groups <- unique(cellInfo$Group)
      if (length(groups) < 2) {
        message("Only one clinical group present; skipping regulon heatmap.")
      } else {
        group_means <- sapply(groups, function(g) {
          rowMeans(auc_mat[sig_regulons, rownames(cellInfo)[cellInfo$Group == g], drop = FALSE])
        })
        colnames(group_means) <- groups
        keep <- apply(group_means, 1, function(r) sd(r) > 0) # drop NaN-prone flat rows
        group_means <- group_means[keep, , drop = FALSE]
        gm_z <- t(scale(t(group_means)))
        rownames(gm_z) <- gsub(" \\(.*\\)", "", rownames(gm_z))
        pdf(file.path(out_dir, "heatmap_sig_regulons.pdf"),
          width = 7,
          height = max(6, nrow(gm_z) * 0.22 + 2)
        )
        pheatmap(gm_z,
          color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
          clustering_method = "ward.D2", fontsize_row = 8, fontsize_col = 10,
          border_color = NA,
          main = sprintf(
            "Top %d significant regulons: mean AUC z-score per group",
            nrow(gm_z)
          )
        )
        dev.off()
      }
    },
    error = function(e) message("regulon heatmap failed (non-fatal): ", conditionMessage(e))
  )
  invisible(res)
}

# Full SCENIC run for a (sub)object, isolated in run_dir, then differential.
run_scenic <- function(obj, spec, run_dir = file.path(SCENIC_RUN_DIR, spec$name),
                       out_dir = file.path(OUTPUT_ROOT, spec$name, "scenic")) {
  suppressPackageStartupMessages({
    library(Seurat)
    library(SCENIC)
    library(AUCell)
    library(GENIE3)
    library(RcisTarget)
  })
  run_dir <- normalizePath(run_dir, mustWork = FALSE)
  db_dir <- normalizePath(SCENIC_DB_DIR, mustWork = FALSE)
  dir.create(file.path(run_dir, "int"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(run_dir, "output"), recursive = TRUE, showWarnings = FALSE)

  exprMat <- as.matrix(GetAssayData(obj, layer = "counts"))
  exprMat <- exprMat[rowSums(exprMat > 0) >= 0.01 * ncol(exprMat), ]
  message(sprintf("SCENIC %s: %d cells, %d genes", spec$name, ncol(exprMat), nrow(exprMat)))

  cellInfo <- data.frame(
    Condition = as.character(obj$condition),
    Group = trimws(as.character(obj$Group)),
    row.names = colnames(exprMat)
  )

  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(run_dir)
  saveRDS(cellInfo, "./int/cellInfo.Rds")

  data(defaultDbNames)
  dbs <- defaultDbNames[["hgnc"]]
  # initializeScenic resolves motifAnnotations_hgnc by name; recent RcisTarget ships it as _v9 (or motifAnnotations), so load whichever exists into globalenv.
  if (!exists("motifAnnotations_hgnc", envir = globalenv())) {
    anno <- tryCatch(
      {
        data("motifAnnotations_hgnc_v9", package = "RcisTarget", envir = environment())
        get("motifAnnotations_hgnc_v9", envir = environment())
      },
      error = function(e) {
        data("motifAnnotations", package = "RcisTarget", envir = environment())
        get("motifAnnotations", envir = environment())
      }
    )
    assign("motifAnnotations_hgnc", anno, envir = globalenv())
  }
  scenicOptions <- initializeScenic(
    org = "hgnc", dbDir = db_dir, dbs = dbs,
    nCores = SCENIC_N_CORES, datasetTitle = spec$name
  )
  scenicOptions@inputDatasetInfo$cellInfo <- "./int/cellInfo.Rds"
  saveRDS(scenicOptions, "./int/scenicOptions.Rds")

  genesKept <- geneFiltering(exprMat, scenicOptions)
  exprMat_filtered <- exprMat[genesKept, ]
  runCorrelation(exprMat_filtered, scenicOptions)
  runGenie3(log2(exprMat_filtered + 1), scenicOptions)

  exprMat_all_log <- log2(exprMat + 1)
  # SCENIC step-3/4 internals resolve scenicOptions by name from globalenv; keep a global copy in sync or runSCENIC_4 aborts.
  assign("scenicOptions", scenicOptions, envir = globalenv())
  scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
  scenicOptions <- runSCENIC_2_createRegulons(scenicOptions)
  assign("scenicOptions", scenicOptions, envir = globalenv())
  scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_all_log)
  assign("scenicOptions", scenicOptions, envir = globalenv())
  # Binarization tail is non-fatal; the differential only needs the step-3 regulon AUC (int/3.4_regulonAUC.Rds).
  step4_ok <- tryCatch(
    {
      scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
      TRUE
    },
    error = function(e) {
      message("runSCENIC_4 tail failed (non-fatal): ", conditionMessage(e))
      FALSE
    }
  )
  # Promote to the _final name only if step 4 completed; otherwise save under a step-3 name, not mislabeled.
  if (step4_ok) {
    saveRDS(scenicOptions, "./int/scenicOptions_final.Rds")
  } else {
    saveRDS(scenicOptions, "./int/scenicOptions_step3.Rds")
    message("Saved pre-binarize object as int/scenicOptions_step3.Rds")
  }
  setwd(old_wd)

  scenic_differential(run_dir, spec, out_dir)
  invisible(run_dir)
}
