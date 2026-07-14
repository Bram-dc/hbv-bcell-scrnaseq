# Transcriptional signatures of ongoing HBV infection in HBcAg-specific B cells in chronic hepatitis B

Analysis code for my bachelor thesis in Nanobiology at Erasmus MC. Both reports are included here:
the [thesis](thesis.pdf) (HBcAg-specific B cells) and a companion
[agentic pipeline-demonstration report](agentic-pipeline-report.pdf) applying the same pipeline to
HBsAg-specific B cells.

## Abstract

Chronic hepatitis B remains a major global health problem, yet B cells are studied far less than
T cells. Hepatitis B virus (HBV)-specific B cells behave dysfunctionally during chronic infection.
Most prior transcriptional analyses profile the total B-cell population rather than the rare cells
that recognize the HBV antigens. This study asks how the transcriptional state of hepatitis B core
antigen (HBcAg)-specific B cells relates to a chronic HBV infection: whether infection leaves a
distinct signature on the transcriptional profile, whether it differs between B-cell subsets, and how
it tracks the viral activity across immune-active (IA), immune-control (IC), and functional-cure (FC)
states. A second aim was to construct a reusable configuration-driven analysis pipeline for these
comparisons. Ongoing infection was expected to leave a signature, strongest in immune-active
patients. Existing plate-based single-cell RNA sequencing of fluorescence-activated cell sorting
(FACS)-sorted cells was analyzed for quality control, clustering, differential expression, gene-set
enrichment, per-cell (single-sample) scoring, and transcription-factor analysis. The results suggest
an activation, interferon, and antigen-presentation program, strongest in IA, weaker in IC, and
lowest in FC, and shared across the subsets rather than confined to one, contrary to the initial
expectation. A preliminary cross-virus comparison with hepatitis C virus (HCV)-specific B cells
indicated that virus identity was not the main source of this variation. Overall, HBcAg-specific B
cells appear to carry a transcriptional signature of activation, interferon, and antigen presentation
that tracks the level of ongoing HBV viral activity.

## The code

Starting from a raw Seurat object of sorted, plate-based scRNA-seq, the code runs QC, clustering,
differential expression, gene-set enrichment (GSEA and ssGSEA), curated module scores, and SCENIC
across a few clinical-state comparisons (IA, IC and FC HBcAg-specific cells, plus an HCV-E2 arm).
It is split three ways:

- `config/` holds the constants (`config.R`) and the comparison registry (`comparisons.R`).
  Edit these to change what gets compared.
- `R/` holds the analysis functions, one file per step (data, subset, cluster, deg, gsea,
  ssgsea, modules, markers, plots, batch, scenic). These are sourced, never run directly.
- `pipeline/` holds the numbered scripts that source the above and write results to `output/`.

## Running

Point `DATA_FILE` in `config/config.R` at the raw Seurat object (default
`data/full_annotation.RData`, object `merged_data`), then from the repo root:

```sh
Rscript pipeline/run_all.R
```

This builds the QC'd, clustered cache once (`00_prepare`), then runs steps `01`-`06` and `08`.
Any step also runs on its own, e.g. `Rscript pipeline/03_deg.R`. SCENIC is heavy and gated off;
run it per comparison with `RUN_SCENIC=TRUE Rscript pipeline/07_scenic.R hbcag_ia_vs_fc`. The
signature score (`09_signature.R`) is a separate post-hoc step.

You need R with Seurat v5 and the usual analysis packages (tidyverse, clusterProfiler,
enrichplot, GSVA, msigdbr, DoubletFinder, org.Hs.eg.db; plus SCENIC, AUCell, GENIE3 and
RcisTarget for step 07). `data/` and `output/` are git-ignored, so the pipeline regenerates
everything under `output/`.

## Acknowledgments

I would like to thank my supervisors for their guidance during this project. I am grateful to Harmen
van de Werken of the Department of Immunology for supervising the work, to Lucas Brock for his close
day-to-day supervision, and to André Boonstra for hosting the project in the Department of
Gastroenterology and Hepatology at Erasmus MC. I am also grateful to Dwin Grashof for his weekly
feedback during the group meetings. Finally, I want to thank the Boonstra laboratory for generating
the single-cell data analyzed here, as well as everyone involved in sample collection and sequencing.

During this thesis, I have used a large language model (Claude Opus, Anthropic) in combination with an
inline AI-driven grammar check (Grammarly) to assist with drafting and editing the text and to give
suggestions for clarity and language. It has also assisted in writing and debugging parts of the
analysis code. It was not used to generate the underlying data or results. It has also been used to
quickly identify genes of interest from large sets of outputs. I reviewed and revised these results
before using them in the analysis.
