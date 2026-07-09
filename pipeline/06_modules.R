# Curated module scores (T-help receptivity, activation/AP-1, ISG, TNF-NF-kB, atypical) per comparison, whole + per cluster. -> output/<cmp>/modules/

source("config/config.R")
source("config/comparisons.R")
invisible(lapply(list.files("R", "[.]R$", full.names = TRUE), source))

obj <- prepare()

for (spec in comparisons_with("modules")) {
  sub <- subset_comparison(obj, spec)
  run_modules(sub, spec)
}
message("\nModule scores complete.")
