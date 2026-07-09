# GSEA from the DEG output tree for comparisons opting into "gsea". -> output/<cmp>/gsea/

source("config/config.R")
source("config/comparisons.R")
invisible(lapply(list.files("R", "[.]R$", full.names = TRUE), source))

term2gene <- term2gene_combined()

for (spec in comparisons_with("gsea")) {
  run_gsea_from_deg(spec, term2gene)
}
message("\nGSEA complete.")
