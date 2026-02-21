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

