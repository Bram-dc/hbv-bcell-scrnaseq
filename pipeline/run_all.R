# Run each pipeline step in a fresh Rscript process, logging to logs/.

source("config/config.R")
source("config/comparisons.R")

dir.create("logs", showWarnings = FALSE)

run_step <- function(script, args = character(0)) {
  tag <- gsub("[^A-Za-z0-9_]", "_", paste(c(script, args), collapse = "_"))
  log <- file.path("logs", paste0(tag, ".log"))
  message(sprintf("\n>>> %s %s", script, paste(args, collapse = " ")))
  status <- system2("Rscript",
    c(file.path("pipeline", paste0(script, ".R")), args),
    stdout = log, stderr = log
  )
  if (status != 0) {
    stop(sprintf("Step failed (status %d): %s - see %s", status, script, log))
  }
  message(sprintf("    done -> %s", log))
}

steps <- c(
  "00_prepare", "01_describe", "02_cluster",
  "03_deg", "04_gsea", "05_ssgsea", "06_modules", "08_synthesis"
)
for (s in steps) run_step(s)

if (isTRUE(RUN_SCENIC)) {
  for (spec in comparisons_with("scenic")) run_step("07_scenic", spec$name)
} else {
  message("\nSCENIC skipped (RUN_SCENIC = FALSE).")
}

message("\n=== run_all complete ===")
