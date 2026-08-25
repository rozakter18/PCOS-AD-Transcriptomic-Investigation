# PCOS–AD Transcriptomic Investigation

A multi-tissue, multi-method transcriptomic investigation of potential molecular links between **Polycystic Ovary Syndrome (PCOS)** and **Alzheimer’s disease (AD)** using publicly available GEO gene-expression datasets.

> **Repository status:** Private research repository  
> **Author:** Roza Akter  
> **Institution:** BRAC University

## Overview

This project investigates whether PCOS gene-expression profiles contain a detectable molecular signal related to Alzheimer’s disease risk.

The analysis began with an individual-level machine-learning approach: train an AD-vs-Control classifier and project PCOS samples into that model to estimate an “AD-likeness” score. After diagnosing cross-platform normalization issues and repeating the analysis in a tissue- and platform-matched blood setting, the classifier itself proved valid but did not identify a reliable PCOS-favoring individual-level AD signal.

The project then shifted to pathway-level analysis using rank-based GSEA, which revealed a more consistent population-level pattern across independent PCOS cohorts.

## Research Questions

1. Can an AD-trained gene-expression classifier detect an individual-level AD-like signal in PCOS samples?
2. Does the result remain stable across different tissues and platforms?
3. Do PCOS and AD share convergent or divergent pathway-level molecular signatures?

## Datasets

| GEO Accession | Disease / Comparison | Tissue | Platform | Samples |
|---|---|---|---|---:|
| **GSE33000** | AD vs. Control | Brain (PFC) | Microarray | 467 |
| **GSE84958** | PCOS vs. Control | Adipose | RNA-seq | 13 used |
| **GSE54248** | PCOS vs. Control | Peripheral blood | Illumina HT-12 v4 | 8 |
| **GSE63061** | AD vs. Control | Peripheral blood | Illumina HT-12 v4 | 273 |

Two technical outliers from GSE84958 were excluded following PCA-based diagnosis.

## Analysis Workflow

### 1. Initial AD Classification and PCOS Projection

- Train a logistic-regression classifier on AD vs. Control expression data
- Reduce dimensionality using PCA
- Project PCOS samples into the same latent space
- Interpret classifier output as an exploratory “AD-likeness” score

### 2. Cross-Platform Harmonization

- Restrict datasets to common genes
- Apply **ComBat** for batch correction
- Re-run the classifier after harmonization

### 3. Honest Model Validation

- Stratified 5-fold cross-validation
- PCA refit within each training fold
- Out-of-fold predictions used for AUC and accuracy
- Avoids information leakage from held-out samples

### 4. Tissue-Matched Replication

A second analysis was performed using blood data on both sides:

- **AD:** GSE63061
- **PCOS:** GSE54248
- Same tissue
- Same microarray platform

This removes major tissue and platform confounding from the original brain-vs-adipose comparison.

### 5. Independent Projection Check

A model-free nearest-centroid method using Pearson correlation was used as an independent check of the classifier-based projection.

### 6. Differential Expression and Outlier Diagnosis

- DESeq2-based analysis for GSE84958
- PCA-based outlier detection
- Cook’s-distance and convergence diagnostics
- Single-gene inference was treated cautiously because of the very small PCOS sample size

### 7. Pathway-Level Analysis

Rank-based GSEA was performed using:

- **MSigDB Hallmark**
- **KEGG Legacy**
- **fgsea**

Genes were ranked using:

- limma moderated t-statistics for microarray datasets
- DESeq2 Wald statistics for RNA-seq

## Key Results

### Individual-Level Classification

The AD classifiers themselves were valid:

- **Brain-based classifier:** CV AUC = **0.963**
- **Blood-based classifier:** CV AUC = **0.755**

However, PCOS samples did **not** show higher AD-likeness than controls in either tissue setting.

| Projection | PCOS Mean | Control Mean |
|---|---:|---:|
| Brain-based classifier | 0.72 | 0.85 |
| Blood-based classifier | 0.43 | 0.56 |
| Blood nearest-centroid | 0.45 | 0.52 |

This indicates that the available public datasets do not support reliable individual-level AD-risk classification from PCOS expression profiles.

### Pathway-Level Findings

A replicated pathway-level signal was detected across independent PCOS cohorts.

Among pathways enriched in both PCOS datasets, the dominant pattern involved **immune and inflammatory pathways being suppressed in PCOS while elevated in AD brain**.

Representative pathways include:

- TNF-α / NF-κB signaling
- IL6 / JAK / STAT3 signaling
- Interferon-γ response
- Inflammatory response
- Complement
- Chemokine signaling
- Leukocyte transendothelial migration

A strong ribosome-related enrichment signal also appeared in both PCOS datasets, but it is treated cautiously because ribosomal pathways can be technically confounded in transcriptomic studies.

## Repository Structure

```text
PCOS_AD/
│
├── GSE33000/
├── GSE54248/
├── GSE63061/
├── GSE84958/
├── results/
│
├── combat_harmonization_classifier.R
├── gsea_ad_pcos.R
├── gsea_three_way_comparison.R
├── PCOS_AD_Investigation_Report.pdf
└── README.md
```

## Main Scripts

### `combat_harmonization_classifier.R`

Handles cross-dataset harmonization, classifier training, validation, and PCOS projection.

### `gsea_ad_pcos.R`

Runs pathway-level enrichment analysis comparing AD and PCOS expression signatures.

### `gsea_three_way_comparison.R`

Compares enrichment patterns across AD brain, PCOS adipose, and PCOS blood datasets.

## Main Findings

- Cross-platform scale mismatch can create misleading classifier outputs if data are projected without harmonization.
- Correctly cross-validated AD classifiers perform well, but this does not imply that PCOS samples contain an individual-level AD-risk signal.
- The null individual-level result replicated across tissue-matched and model-free analyses.
- Pathway-level analysis produced a stronger and more reproducible biological signal than individual-level classification.
- Immune/inflammatory pathways showed a replicated directionally divergent pattern between PCOS and AD.

## Limitations

- Very small PCOS cohort sizes
- Cross-sectional rather than longitudinal data
- Strong assumptions required for cross-study ComBat harmonization
- Considerable gene loss during cross-platform matching
- Only one AD dataset per tissue
- No direct longitudinal PCOS-to-AD outcome data

## Future Work

Potential extensions include:

- pooling additional PCOS cohorts
- replication in additional tissue types
- leading-edge gene analysis of significant pathways
- independent AD-side replication
- longitudinal cohorts linking PCOS diagnosis to later AD outcomes

## Report

The complete investigation is documented in:

[`PCOS_AD_Investigation_Report.pdf`](PCOS_AD_Investigation_Report.pdf)

## Data Availability

All four datasets used in this project are publicly available from the **NCBI Gene Expression Omnibus (GEO)**:

- GSE33000
- GSE84958
- GSE54248
- GSE63061

Raw and processed data included locally should be reviewed before pushing to GitHub, especially for large files.

## Academic Note

This repository documents an exploratory computational-biology investigation. The individual-level classifier result is negative, while the pathway-level results are hypothesis-generating and should not be interpreted as clinical evidence of Alzheimer’s risk in individual PCOS patients.
