# Whole-dataset description of cell and patient composition by group and antigen, plus batch diagnostics. -> output/describe/

source("config/config.R")
source("config/comparisons.R")
invisible(lapply(list.files("R", "[.]R$", full.names = TRUE), source))

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

out_dir <- file.path(OUTPUT_ROOT, "describe")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

obj <- prepare()
meta <- obj@meta.data |> mutate(Group = trimws(Group), Antigen = trimws(Antigen))

group_colours <- c(
  "IA" = "#E63946", "SNEG Ab+" = "#457B9D", "SNEG Ab-" = "#A8DADC",
  "Acute" = "#F4A261", "HCV" = "#2A9D8F", "ENEG" = "#264653", "HC" = "#E9C46A",
  "NUC" = "#8338EC", "SARS" = "#FB5607", "Memory" = "#AECBFA", "Naive" = "#D4E09B",
  "IC(high)" = "#C77DFF", "IC(low)" = "#E0AAFF", "IT" = "#606C38",
  "Atypical memory" = "#BC6C25", "Active memory" = "#DDA15E",
  "Resting memory" = "#FEFAE0", "T-cell" = "#B5838D", "Total B cell" = "#6D6875"
)
antigen_colours <- c(
  "HBcAg" = "#E63946", "HBsAg" = "#457B9D", "E2" = "#2A9D8F", "SARS-COV2" = "#F4A261",
  "Naive" = "#D4E09B", "Memory" = "#AECBFA", "Total B cell" = "#6D6875",
  "Atypical memory" = "#BC6C25", "Active memory" = "#DDA15E",
  "Resting memory" = "#FEFAE0", "T-cell" = "#B5838D"
)

save_plot <- function(p, name, w, h) {
  ggsave(file.path(out_dir, paste0(name, ".png")), p,
    width = w, height = h,
    dpi = 200, bg = "white"
  )
  message("Saved: ", name)
}

save_plot(
  ggplot(count(meta, Group), aes(reorder(Group, n), n, fill = Group)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = n), hjust = -0.15, size = 3) +
    scale_fill_manual(values = group_colours, na.value = "grey70") +
    coord_flip() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(title = "Cells per Group", x = NULL, y = "Cell count") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none"),
  "cells_per_group", 8, 6
)

save_plot(
  ggplot(count(meta, Antigen), aes(reorder(Antigen, n), n, fill = Antigen)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = n), hjust = -0.15, size = 3) +
    scale_fill_manual(values = antigen_colours, na.value = "grey70") +
    coord_flip() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(title = "Cells per Antigen", x = NULL, y = "Cell count") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none"),
  "cells_per_antigen", 8, 6
)

cross <- meta |>
  count(Group, Antigen) |>
  complete(Group, Antigen, fill = list(n = 0))
save_plot(
  ggplot(cross, aes(Antigen, Group, fill = log1p(n))) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = ifelse(n > 0, n, "")), size = 2.8) +
    scale_fill_gradient(low = "white", high = "#E63946", name = "log1p(n)") +
    labs(title = "Group versus Antigen composition", x = "Antigen", y = "Group") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1)),
  "group_x_antigen", 10, 7
)

prop <- meta |>
  count(Group, Antigen) |>
  group_by(Group) |>
  mutate(prop = n / sum(n))
save_plot(
  ggplot(prop, aes(Group, prop, fill = Antigen)) +
    geom_col(position = "stack", width = 0.75) +
    scale_fill_manual(values = antigen_colours, na.value = "grey70") +
    scale_y_continuous(labels = scales::percent) +
    labs(title = "Antigen composition within each Group", x = NULL, y = "Proportion") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1), legend.key.size = unit(0.4, "cm")),
  "antigen_proportion", 10, 6
)

sample_col <- if ("Sample" %in% colnames(meta)) "Sample" else "orig.ident"
pm <- meta |> rename(Patient = !!sym(sample_col))

save_plot(
  ggplot(
    count(pm, Patient, Group) |> group_by(Patient) |> mutate(total = sum(n)),
    aes(reorder(Patient, total), n, fill = Group)
  ) +
    geom_col(width = 0.75) +
    scale_fill_manual(values = group_colours, na.value = "grey70") +
    coord_flip() +
    labs(
      title = "Cells captured per patient", subtitle = "Coloured by Group",
      x = NULL, y = "Cell count"
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.key.size = unit(0.4, "cm")),
  "patient_cells_by_group", 10, 8
)

save_plot(
  ggplot(
    pm |> distinct(Patient, Group) |> count(Group),
    aes(reorder(Group, n), n, fill = Group)
  ) +
    geom_col(width = 0.7) +
    geom_text(aes(label = n), hjust = -0.2, size = 3.5) +
    scale_fill_manual(values = group_colours, na.value = "grey70") +
    coord_flip() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
    labs(title = "Patients per Group", x = NULL, y = "Number of patients") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none"),
  "patients_per_group", 8, 6
)

save_plot(
  ggplot(count(pm, Patient, Group), aes(Group, n, fill = Group)) +
    geom_boxplot(outlier.shape = 21, width = 0.55, alpha = 0.8) +
    geom_jitter(width = 0.15, size = 1.5, alpha = 0.6) +
    scale_fill_manual(values = group_colours, na.value = "grey70") +
    labs(title = "Cell count distribution per patient within Group", x = NULL, y = "Cells per patient") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1), legend.position = "none"),
  "cellcount_distribution", 10, 6
)

batch_diagnostics(obj, out_dir)
message("\nDescribe outputs saved to: ", out_dir)
