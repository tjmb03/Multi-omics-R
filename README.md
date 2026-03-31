> **Disclaimer:** All datasets in this repository are simulated or 
> pseudodata generated solely for methodological demonstration purposes. 
> No proprietary, confidential, patient-derived, or employer-affiliated 
> data is included. This work represents independent research and 
> educational development conducted outside of any employment context 
> and does not reflect the proprietary methods, data, or intellectual 
> property of any employer or collaborator.
> 
> This repository is released under the [MIT License](LICENSE).
> © 2026 Bo Ma (tjmb03). Reuse with attribution.

# Multi-Omics Integration Toolkit (R)

A modular framework for multi-omics integration, biomarker discovery, and disease stratification.

This repository implements structured workflows for:

- Latent factor modeling (MOFA2)
- Supervised multi-block integration (DIABLO / mixOmics)
- Joint clustering (iCluster)
- Network-based module discovery (WGCNA)
- Bayesian network inference
- Survival association modeling

Designed as a professional systems biology toolkit rather than standalone scripts.

---

# Multi-Omics Analysis Suite — Shiny App

A unified, dark-themed Shiny web application for running five multi-omics integration workflows on paired transcriptomics and proteomics data, with live progress tracking, parameter controls, and HTML/PDF report export.

See its [folder](https://github.com/tjmb03/Multi-omics-R/tree/main/Multi-Omics%20Analysis%20Suite%20%E2%80%94%20ShinyApp) for scripts and usage details.
---

# Core Modules

## 1️⃣ Latent Factor Modeling (MOFA2)

Implements unsupervised multi-view factor analysis to:

- Decompose cross-omics variance
- Identify shared and modality-specific signals
- Prioritize biomarker candidates
- Discover disease-associated latent structure

Applications:
- Mechanism-driven biomarker discovery
- Patient stratification
- Cross-platform signal extraction

---

## 2️⃣ Supervised Multi-Block Integration (DIABLO)

Implements mixOmics DIABLO for:

- Cross-omics correlated feature selection
- Supervised classification
- Integrated biomarker panel discovery
- Multi-block predictive modeling

Applications:
- Translational signature development
- Disease subtype discrimination
- Multi-omic panel construction

---

## 3️⃣ Joint Clustering (iCluster)

Implements integrative clustering to:

- Identify molecular subtypes
- Jointly model multiple omics layers
- Discover cross-omics latent structure

Applications:
- Disease stratification
- Cohort segmentation
- Molecular subtype discovery

---

## 4️⃣ Network-Based Module Discovery (WGCNA)

Implements weighted co-expression network analysis to:

- Identify co-regulated gene modules
- Relate modules to clinical traits
- Associate modules with survival outcomes

Applications:
- Systems-level biomarker prioritization
- Prognostic module discovery
- Network-level interpretation

---

## 5️⃣ Bayesian Network Inference

Implements probabilistic graphical modeling to:

- Infer conditional dependencies
- Explore regulatory relationships
- Generate mechanistic hypotheses

Applications:
- Causal structure exploration
- Regulatory modeling
- Systems biology hypothesis generation

---

# Capabilities Demonstrated

- Multi-view latent factor modeling
- Supervised multi-block integration
- Integrative clustering
- Network analysis
- Survival modeling
- Probabilistic graphical modeling
- TCGA-scale dataset handling
- Translational biomarker strategy

---

# Technical Stack

- R ≥ 4.2
- MOFA2
- mixOmics
- iClusterPlus
- WGCNA
- bnlearn
- survival
- tidyverse

Install core packages:

```r
install.packages(c(
  "mixOmics",
  "WGCNA",
  "bnlearn",
  "survival",
  "tidyverse"
))
```
> © 2026 tjmb03. This project is provided for educational and methodological
demonstration purposes. Source code for the interactive dashboards is **available on request** for academic and research use.
