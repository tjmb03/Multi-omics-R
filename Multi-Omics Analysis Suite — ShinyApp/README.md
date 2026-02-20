# Multi-Omics Analysis Suite — Shiny App

A unified, dark-themed Shiny web application for running five multi-omics integration workflows on paired transcriptomics and proteomics data, with live progress tracking, parameter controls, and HTML/PDF report export.


[![Launch App](https://img.shields.io/badge/Launch-Live%20Shiny%20App-blue?style=for-the-badge)](https://tjmb03.shinyapps.io/MultiOmics_Shiny_App/)


## Overview

The app provides a single browser-based interface to run five published multi-omics methods on a glaucoma dataset (24 samples, transcriptomics + proteomics). All heavy computation runs in a background R process (`callr::r_bg`) so the UI stays fully responsive during analysis. Progress is streamed live to the browser every second.

**Key features:**
- Five workflow selector buttons with live parameter grids
- File upload slots for all five CSV inputs (up to 100 MB each)
- Live progress bar with named stages and percentage display
- Console log streaming render output in real time
- HTML report download (code hidden, results only)
- Optional PDF report download (requires LaTeX/tinytex)
- Generated R script download for manual reproduction in RStudio

---

## Workflows

| # | Workflow | Type | Method |
|---|---|---|---|
| 1 | **MOFA2** | Unsupervised · Dimensionality Reduction | Multi-Omics Factor Analysis — learns shared latent factors across RNA-seq and proteomics |
| 2 | **DIABLO** | Supervised · Biomarker Discovery | Data Integration Analysis for Biomarker discovery using Latent cOmponents — discriminative multi-omic signatures |
| 3 | **iClusterPlus** | Unsupervised · Clustering | Integrative clustering via shared latent variable model with BIC-based cluster selection |
| 4 | **WGCNA** | Unsupervised · Network Analysis | Weighted Gene Co-expression Network Analysis — module detection and hub gene identification |
| 5 | **Bayesian Network** | Probabilistic · Causal Inference | Hill-climbing DAG search with bootstrap edge strength for gene-protein-clinical causal graphs |


## Input Data Files

| File | Description | Required by |
|---|---|---|
| `raw_counts.csv` | RNA-seq raw count matrix (genes x samples) | All workflows |
| `raw_proteomics.csv` | Proteomics abundance matrix (proteins x samples) | MOFA2, DIABLO, iClusterPlus, BayesNet |
| `sample_metadata.csv` | Sample metadata with `sample_id` and diagnosis columns | All workflows |
| `map_transcriptomics.csv` | Gene ID mapping table (Ensembl to symbol) | All workflows |
| `map_proteomics.csv` | Protein ID mapping table | MOFA2, DIABLO, iClusterPlus, BayesNet |

**WGCNA** only requires 3 files: `raw_counts.csv`, `sample_metadata.csv`, `map_transcriptomics.csv`.

---

## Required R Packages

### Shiny App Core

```r
install.packages(c(
  "shiny",      # Web application framework
  "bslib",      # Bootstrap 5 theming + Google Fonts
  "shinyjs",    # JS helpers (enable/disable, runjs, html)
  "callr",      # Background R process for non-blocking render
  "rmarkdown",  # Rmd rendering to HTML and PDF
  "knitr"       # Chunk execution engine
))
```

### Workflow Packages — All in One

```r
# CRAN
install.packages(c(
  "tidyverse",   # Data wrangling (dplyr, tidyr, ggplot2, purrr, etc.)
  "ggplot2",     # Plotting
  "ggrepel",     # Non-overlapping text labels on plots
  "matrixStats", # rowVars, colMedians, etc.
  "pheatmap",    # Heatmap visualisation (DIABLO, iClusterPlus)
  "WGCNA",       # Weighted Gene Co-expression Network Analysis
  "bnlearn",     # Bayesian network structure learning
  "arules",      # Discretisation for Bayesian Network
  "igraph"       # Network graph operations (WGCNA, BayesNet)
))

# Bioconductor
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c(
  "MOFA2",           # Multi-Omics Factor Analysis v2
  "mixOmics",        # DIABLO method
  "iClusterPlus",    # Integrative clustering
  "edgeR",           # TMM normalisation of raw counts
  "clusterProfiler", # GO and KEGG enrichment analysis
  "org.Hs.eg.db",    # Human gene annotation database
  "Rgraphviz"        # DAG visualisation for Bayesian Network
))
```

## Parameter Reference

### MOFA2

| Parameter | Range | Default | Notes |
|---|---|---|---|
| Number of Latent Factors | 2–20 | 10 | More factors = finer variation, higher overfitting risk |
| Max Training Iterations | 100–2000 | 1000 | Uses ELBO convergence, usually stops early |
| Random Seed | any integer | 42 | Fix for reproducibility |

### DIABLO

| Parameter | Range | Default | Notes |
|---|---|---|---|
| Top Variable Genes | 100–5000 | 1000 | Pre-filter before tuning |
| Number of Components | 1–5 | 2 | 2 components sufficient for binary outcome |
| Design Matrix Off-diagonal | 0.1–1.0 | 0.1 | 0.1 prioritises discrimination |
| CV Folds | 3–10 | 5 | Use 3–5 folds with n=24 |
| CV Repeats | 1–20 | 5 | More repeats reduce error variance |

### iClusterPlus

| Parameter | Range | Default | Notes |
|---|---|---|---|
| Max Clusters (k) | 2–5 | 3 | k > 3 unreliable with n=24 |
| Top Variable Genes | 100–1000 | 500 | Genes retained for clustering |
| Lambda Grid Points | 8, 13, 21, 34, 55 | 8 | Must be Fibonacci numbers |
| Max EM Iterations | 10–100 | 20 | 20 sufficient for exploration |

### WGCNA

| Parameter | Range | Default | Notes |
|---|---|---|---|
| Minimum Module Size | 5–100 | 10 | Keep 10–20 with small n |
| Merge Cut Height | 0.10–0.50 | 0.25 | Lower = more merging |
| Deep Split | 0–4 | 2 | Higher = more aggressive splitting |
| Soft Power (0 = auto) | 0–20 | 0 | Auto-selects power for scale-free R2 >= 0.80 |

### Bayesian Network

| Parameter | Range | Default | Notes |
|---|---|---|---|
| Top Variable Genes | 5–50 | 20 | BN runtime is O(2^n) — keep <= 20 |
| Top Variable Proteins | 5–30 | 10 | Combined nodes should not exceed 30 |
| Discretisation Bins | 2–5 | 3 | Fewer bins = more stable with small n |
| Bootstrap Resamples | 10–500 | 50 | R=50 finishes in ~1–2 min with BIC |
| Structure Score | BIC / BDe | BIC | BIC is ~10x faster |
| HC Restarts | 1–10 | 3 | More restarts reduce local-optimum risk |


## Output and Downloads

| Output | Format | Code | Description |
|---|---|---|---|
| HTML Report | `.html` | Hidden | Self-contained report with all plots and tables |
| PDF Report | `.pdf` | Hidden | Same content (requires tinytex) |
| R Script | `.R` | — | Parameter injection script for manual knitting |

Reports are named `WORKFLOW_report_YYYYMMDD_HHMMSS.html/.pdf`.

## Tested With

- R 4.3+, RStudio 2023+
- shiny 1.8+, bslib 0.7+, shinyjs 2.1+, callr 3.7+
- rmarkdown 2.26+, knitr 1.45+
- MOFA2 1.12+ (Bioconductor 3.18)
- mixOmics 6.26+, iClusterPlus 1.38+
- WGCNA 1.72+, bnlearn 4.9+

