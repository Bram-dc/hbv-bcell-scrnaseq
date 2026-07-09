# Build the QC'd, doublet-filtered, clustered cache once for later steps to load.

source("config/config.R")
source("config/comparisons.R")
invisible(lapply(list.files("R", "[.]R$", full.names = TRUE), source))

obj <- prepare(apply_qc = TRUE, use_cache = TRUE)
message(sprintf("Prepared object: %d cells, %d genes", ncol(obj), nrow(obj)))
message(sprintf("HBcAg cells: %d", sum(trimws(obj$Antigen) == "HBcAg")))
