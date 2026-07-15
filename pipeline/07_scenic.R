# SCENIC per comparison/cluster ("diff" = differential only, "loom" = SCope loom). -> output/<cmp>/scenic[/<cluster>]/

source("config/config.R")
source("config/comparisons.R")
invisible(lapply(list.files("R", "[.]R$", full.names = TRUE), source))

args <- commandArgs(trailingOnly = TRUE)
cmp <- if (length(args) >= 1) args[1] else "hbcag_ia_vs_fc"
diff_only <- "diff" %in% args
loom_only <- "loom" %in% args
rest <- setdiff(args, c(cmp, "diff", "loom", "force"))
cl <- if (length(rest) >= 1) rest[1] else NA

spec <- COMPARISONS[[cmp]]
if (is.null(spec)) stop("Unknown comparison: ", cmp)

suffix <- if (is.na(cl)) "" else paste0("_", cl)
run_dir <- file.path(SCENIC_RUN_DIR, paste0(cmp, suffix))
out_dir <- file.path(OUTPUT_ROOT, cmp, "scenic")
if (!is.na(cl)) out_dir <- file.path(out_dir, cl)

if (diff_only) {
  message("SCENIC differential-only: ", cmp, suffix, " (run_dir = ", run_dir, ")")
  scenic_differential(run_dir, spec, out_dir)
  message("Differential written to ", out_dir)
  quit(save = "no", status = 0)
}

# Reads the cached run, so no RUN_SCENIC gate or cisTarget databases needed.
if (loom_only) {
  message("SCENIC loom export: ", cmp, suffix, " (run_dir = ", run_dir, ")")
  scenic_export_loom(run_dir, spec, prepare(), overwrite = "force" %in% args)
  quit(save = "no", status = 0)
}

if (!isTRUE(RUN_SCENIC)) {
  message("RUN_SCENIC is FALSE in config/config.R - skipping SCENIC for ", cmp)
  quit(save = "no", status = 0)
}

obj <- prepare()
sub <- subset_comparison(obj, spec)
if (!is.na(cl)) sub <- with_clusters(sub, cl)

run_scenic(sub, spec, run_dir = run_dir, out_dir = out_dir)
message("\nSCENIC complete: ", cmp, suffix)
