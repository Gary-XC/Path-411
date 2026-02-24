# Group Assignment 1 – Genomic Data Preprocessing Pipeline

## Objective

Build a full preprocessing workflow for publicly available miRNA-seq lymphoma data.

---

## Dataset

ICGC MALY-DE miRNA dataset  
Filtered to Diffuse Large B-Cell Lymphoma (DLBCL)

---

## Pipeline Steps

- Data ingestion and restructuring
- Sample filtering based on diagnosis codes
- Relative frequency normalization
- Log2 transformation
- Outlier detection (IQR + Spearman correlation)
- Batch effect inspection
- Quantile-based low expression filtering
- Final QC visualization

---

##  Visualizations Generated

- Total counts per sample
- IQR distributions
- Correlation-based outlier detection plots
- Post-normalization comparisons

---

##  ML Relevance

This mirrors real-world preprocessing in ML:

- Data cleaning
- Feature filtering
- Outlier removal
- High-dimensional dataset conditioning