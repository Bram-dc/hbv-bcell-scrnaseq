# Gene-set and gene-list helpers for HGNC, protein-coding, and MSigDB sets.

suppressPackageStartupMessages({
  library(msigdbr)
})

# Load and cache the HGNC complete gene set.
load_hgnc <- function(cache = HGNC_CACHE, url = HGNC_URL) {
  if (file.exists(cache)) {
    return(readRDS(cache))
  }
  message("Downloading HGNC complete gene set...")
  hgnc <- read.delim(url, check.names = FALSE)
  dir.create(dirname(cache), recursive = TRUE, showWarnings = FALSE)
  saveRDS(hgnc, cache)
  hgnc
}

# Immunoglobulin gene symbols (HGNC locus_type).
ig_genes <- function(hgnc = load_hgnc()) {
  hgnc$symbol[hgnc$locus_type %in%
    c("immunoglobulin gene", "immunoglobulin pseudogene")]
}

# Long non-coding RNA symbols (HGNC locus_type).
lnc_genes <- function(hgnc = load_hgnc()) {
  hgnc$symbol[hgnc$locus_type == "RNA, long non-coding"]
}

# Protein-coding subset of symbols (org.Hs.eg.db GENETYPE).
protein_coding_genes <- function(symbols) {
  gt <- suppressMessages(AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = unique(symbols), columns = "GENETYPE", keytype = "SYMBOL"
  ))
  unique(gt$SYMBOL[!is.na(gt$GENETYPE) & gt$GENETYPE == "protein-coding"])
}

# Load one MSigDB collection.
load_gene_sets <- function(collection, subcategory = NULL) {
  sets <- msigdbr(
    species = "Homo sapiens",
    collection = collection,
    subcollection = subcategory
  )
  term2gene <- sets[, c("gs_name", "gene_symbol")]
  gene_sets <- lapply(split(sets$gene_symbol, sets$gs_name), unique)

  label <- dplyr::case_when(
    collection == "H" ~ "Hallmark",
    collection == "C2" & !is.null(subcategory) ~ paste0("C2 ", subcategory),
    collection == "C2" ~ "C2 Canonical",
    collection == "C7" & !is.null(subcategory) ~ paste0("C7 ", subcategory),
    collection == "C7" ~ "C7 Immunologic",
    TRUE ~ paste(collection, subcategory)
  )
  slug <- tolower(gsub("[^a-zA-Z0-9]+", "_", trimws(label)))
  message(sprintf("Loaded %d gene sets  [%s]", length(gene_sets), label))

  list(
    term2gene = term2gene, gene_sets = gene_sets,
    label = label, slug = slug, n_sets = length(gene_sets)
  )
}

# Collections flagged for GSEA and/or ssGSEA.
COLLECTIONS <- list(
  hallmark = list(
    collection = "H", subcategory = NULL,
    use_gsea = TRUE, use_ssgsea = TRUE
  ),
  c2_reactome = list(
    collection = "C2", subcategory = "CP:REACTOME",
    use_gsea = TRUE, use_ssgsea = TRUE
  ),
  c7_immune = list(
    collection = "C7", subcategory = "IMMUNESIGDB",
    use_gsea = TRUE, use_ssgsea = FALSE
  ) # too slow per-cell
)

# Load all collections enabled for the given flag.
get_collections <- function(type = c("gsea", "ssgsea")) {
  type <- match.arg(type)
  flag <- paste0("use_", type)
  keep <- Filter(function(x) isTRUE(x[[flag]]), COLLECTIONS)
  lapply(keep, function(x) load_gene_sets(x$collection, x$subcategory))
}

# Combined term2gene: C7 immune + C2 Reactome + C5 GO:BP.
term2gene_combined <- function() {
  c7 <- msigdbr(species = "Homo sapiens", collection = "C7", subcollection = "IMMUNESIGDB")
  c2 <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME")
  go <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP")
  combined <- dplyr::bind_rows(
    dplyr::select(c7, "gs_name", "gene_symbol"),
    dplyr::select(c2, "gs_name", "gene_symbol"),
    dplyr::select(go, "gs_name", "gene_symbol")
  )
  combined <- as.data.frame(combined)
  message(sprintf("Combined gene sets: %d terms", dplyr::n_distinct(combined$gs_name)))
  combined
}
