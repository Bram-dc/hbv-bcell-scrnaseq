# Score the findings-derived viral-activity signature on the full HBcAg population across states and subsets. -> output/signature/

source("config/config.R")
source("config/comparisons.R")
invisible(lapply(list.files("R", "[.]R$", full.names = TRUE), source))

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

out_dir <- file.path(OUTPUT_ROOT, "signature")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Full HBcAg population, clustered with the canonical recipe.
obj <- prepare()
ag <- cluster_umap(subset_antigen(obj, "HBcAg"))

# Map clinical phases onto the three viral activity states.
g <- trimws(as.character(ag@meta.data[["Group"]]))
state <- ifelse(g == "IA", "IA",
  ifelse(g %in% c("IC(high)", "IC(low)"), "IC",
    ifelse(g %in% FC_GROUPS, "FC", NA_character_)
  )
)
keep <- !is.na(state) & as.character(ag$seurat_clusters) %in% c("0", "1", "2")
state_keep <- state[keep]
ag <- subset(ag, cells = colnames(ag)[keep])
ag$state <- factor(state_keep, levels = c("IA", "IC", "FC"))

# Score the findings-derived signature (genes present only).
genes <- intersect(SIGNATURE_SET, rownames(ag))
message(sprintf("Signature genes present: %d / %d", length(genes), length(SIGNATURE_SET)))
ag <- AddModuleScore(ag, features = list(genes), name = "sig_", seed = 42, search = FALSE)
ag$score <- ag$sig_1
ag$subset <- factor(
  cluster_name(as.character(ag$seurat_clusters), "HBcAg"),
  levels = cluster_name(c("0", "1", "2"), "HBcAg")
)

# Long data: an "All subsets" facet plus one per cluster.
df_all <- data.frame(state = ag$state, score = ag$score, facet = "All subsets")
df_cl <- data.frame(state = ag$state, score = ag$score, facet = as.character(ag$subset))
df <- rbind(df_all, df_cl)
df$facet <- factor(df$facet, levels = c("All subsets", cluster_name(c("0", "1", "2"), "HBcAg")))

# Median score per state and facet.
# Descriptive only; the per-cell test is pseudoreplicated.
summ <- df |>
  group_by(facet, state) |>
  summarise(n_cells = dplyr::n(), median_score = median(score, na.rm = TRUE), .groups = "drop")
write.csv(summ, file.path(out_dir, "signature_score_summary.csv"), row.names = FALSE)
cat("\n=== Median viral-activity-signature score per state and subset ===\n")
print(as.data.frame(summ))

pal <- c("IA" = "#E63946", "IC" = "#F4A261", "FC" = "#457B9D")
p <- ggplot(df, aes(state, score, fill = state)) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.7) +
  geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA, colour = "grey30") +
  facet_wrap(~facet, nrow = 1) +
  scale_fill_manual(values = pal) +
  labs(
    title = "Viral-activity signature score across states and subsets",
    subtitle = sprintf(
      "Findings-derived module (%d genes: interferon, MHC class II, activation) scored per cell; box = IQR/median",
      length(genes)
    ),
    x = NULL, y = "Signature score"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )
ggsave(file.path(out_dir, "signature_by_state_and_subset.png"), p,
  width = 13, height = 4.6, dpi = 200, bg = "white"
)
message(sprintf("\nSignature figure saved to: %s", out_dir))
