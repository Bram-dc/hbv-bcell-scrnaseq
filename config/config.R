# Shared analysis constants: paths, QC, embedding/clustering, method cutoffs, marker/module sets, and cluster labels.

DATA_FILE <- "./data/full_annotation.RData"
OBJECT_NAME <- "merged_data"

CACHE_DIR <- "./output/cache"
SINGLET_CACHE <- file.path(CACHE_DIR, "seurat_singlets.rds")
PREQC_CACHE <- file.path(CACHE_DIR, "seurat_singlets_preqc.rds")
HGNC_CACHE <- file.path(CACHE_DIR, "hgnc_complete_set.rds")
HGNC_URL <- "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"

OUTPUT_ROOT <- "./output"

BATCH_VAR <- "Sequencing_batch"

PCA_DIMS <- 1:20
CLUSTER_RES <- 0.4

QC_MIN_FEATURES <- 200
QC_MAX_FEATURES <- 3000
QC_MAX_MT <- 20
QC_MIN_COUNTS <- 500

DF_PK <- 0.03
DF_PN <- 0.25
DF_NEXP_FRAC <- 0.075

DEG_LOGFC <- 0.25
DEG_MIN_PCT <- 0.10
DEG_P <- 0.05

GSEA_MIN_GS <- 10
GSEA_MAX_GS <- 500
GSEA_PCUT <- 0.25

TOP_N <- 20
TOP_DE_N <- 8
MIN_PCT_EXPRESS <- 10

RUN_SCENIC <- tolower(Sys.getenv("RUN_SCENIC", "FALSE")) %in% c("true", "1", "yes")
SCENIC_DB_DIR <- "./cisTarget_databases"
SCENIC_N_CORES <- as.integer(Sys.getenv("SCENIC_N_CORES", "16"))
SCENIC_RUN_DIR <- "./output/scenic"

VIOLIN_MARKERS <- c(
  "CR2", "CD27", "CD69", "CD83", "FCRL5", "TBX21"
)

MODULE_SETS <- list(
  resting = c(
    "CD27", "CR2", "CXCR5", "TNFRSF13B", "SELL"
  ),
  activation = c(
    "IRF1", "STAT1", "STAT3", "CD69", "CD83", "CD86", "TFRC", "FAS",
    "FOS", "FOSB", "JUN", "JUNB", "JUND"
  ),
  atypical = c(
    "FCRL5", "FCRL3", "ITGAX", "TBX21", "ZEB2", "FCGR2B", "SIGLEC6", "CD72"
  ),
  t_help_receptivity = c(
    "CD40", "CD27", "IL4R", "TNFRSF13B", "TNFRSF13C", "FAS",
    "CD80", "CD86", "HLA-DRA", "CXCR5", "CXCR4", "SELL"
  ),
  isg = c(
    "ISG15", "MX1", "MX2", "IFI44", "IFI44L",
    "OAS1", "LY6E", "IFITM3", "IRF7"
  ),
  tnf_nfkb = c(
    "NFKB1", "NFKBIA", "TNF", "TNFAIP3",
    "TNFRSF1B", "CD83", "FOS", "JUN"
  )
)

ALL_MARKERS <- unique(c(
  "CD19", "MS4A1",
  MODULE_SETS$resting, MODULE_SETS$activation, MODULE_SETS$atypical
))

# Findings-derived viral-activity signature (activation/AP-1, ISG, MHC-II), scored as one module (pipeline/09_signature.R).
SIGNATURE_SET <- unique(c(
  MODULE_SETS$activation, "EGR1", "NR4A1", "NR4A2", "DUSP1",
  MODULE_SETS$isg, "IFI30", "IFI6", "IFITM1",
  "HLA-DRA", "HLA-DRB1", "HLA-DRB5", "HLA-DQA1", "HLA-DQA2", "HLA-DQB1",
  "HLA-DPA1", "HLA-DPB1", "HLA-DMA", "HLA-DMB", "CD74"
))

CLUSTER_LABELS <- list(
  HBcAg = c(
    "0" = "Naive / resting memory",
    "1" = "Activated",
    "2" = "FCRL5+ atypical",
    "3" = "MT-high",
    "4" = "Non-B contaminant"
  )
)

# Map cluster id(s) to their subset name for an antigen, falling back to the id when unlabeled.
cluster_name <- function(id, antigen = "HBcAg") {
  ids <- as.character(id)
  labs <- CLUSTER_LABELS[[antigen[1]]]
  if (is.null(labs)) {
    return(ids)
  }
  ifelse(ids %in% names(labs), unname(labs[ids]), ids)
}
