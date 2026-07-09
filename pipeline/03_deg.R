# Group-placement UMAP for every comparison plus DEG for opt-in comparisons. -> output/<cmp>/umap_condition.png, output/<cmp>/deg/

source("config/config.R")
source("config/comparisons.R")
invisible(lapply(list.files("R", "[.]R$", full.names = TRUE), source))

obj <- prepare()

# For each comparison: resolve the subset, save the group-placement UMAP, then run DEG if it opts in.
for (spec in COMPARISONS) {
  sub <- subset_comparison(obj, spec)
  umap_condition_plot(sub, spec, file.path(OUTPUT_ROOT, spec$name, "umap_condition.png"))
  if ("deg" %in% spec$methods) run_deg(sub, spec)
}
message("\nDEG complete.")
