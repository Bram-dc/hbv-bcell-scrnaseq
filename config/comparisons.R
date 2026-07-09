# The COMPARISONS registry: each entry defines one contrast.

FC_GROUPS <- c("SNEG Ab+", "SNEG Ab-")

ALL_METHODS <- c("deg", "gsea", "ssgsea", "modules", "scenic")

COMPARISONS <- list(
  hbcag_ia_vs_fc = list(
    name = "hbcag_ia_vs_fc",
    antigen = "HBcAg",
    label = "IA vs FC (HBcAg)",
    group_col = "Group",
    ident1 = "IA", groups1 = "IA",
    ident2 = "FC", groups2 = FC_GROUPS,
    per_cluster = TRUE,
    exclude_clusters = c("3", "4"),
    methods = ALL_METHODS
  ),
  hbcag_ic_vs_fc = list(
    name = "hbcag_ic_vs_fc",
    antigen = "HBcAg",
    label = "IC vs FC (HBcAg)",
    group_col = "Group",
    ident1 = "IC", groups1 = c("IC(high)", "IC(low)"),
    ident2 = "FC", groups2 = FC_GROUPS,
    per_cluster = TRUE,
    exclude_clusters = c("3", "4"),
    methods = ALL_METHODS
  ),
  hbcag_ia_vs_ic = list(
    name = "hbcag_ia_vs_ic",
    antigen = "HBcAg",
    label = "IA vs IC (HBcAg)",
    group_col = "Group",
    ident1 = "IA", groups1 = "IA",
    ident2 = "IC", groups2 = c("IC(high)", "IC(low)"),
    per_cluster = TRUE,
    exclude_clusters = c("3", "4"),
    methods = ALL_METHODS
  ),
  hbcag_ia_vs_hcv = list(
    name = "hbcag_ia_vs_hcv",
    antigen = c("HBcAg", "E2"),
    label = "HBcAg IA vs HCV E2",
    group_col = "Group",
    ident1 = "IA", groups1 = "IA",
    ident2 = "HCV", groups2 = "HCV",
    per_cluster = FALSE,
    exclude_clusters = character(0),
    methods = ALL_METHODS
  )
)

comparisons_with <- function(method) {
  Filter(function(s) method %in% s$methods, COMPARISONS)
}
