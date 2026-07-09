# ssGSEA per comparison that opts into "ssgsea". -> output/<cmp>/ssgsea/

source("config/config.R")
source("config/comparisons.R")
invisible(lapply(list.files("R", "[.]R$", full.names = TRUE), source))

obj <- prepare()
collections <- get_collections("ssgsea")

for (spec in comparisons_with("ssgsea")) {
  sub <- subset_comparison(obj, spec)
  run_ssgsea(sub, spec, collections)
}
message("\nssGSEA complete.")
