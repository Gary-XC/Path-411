# Group Assignment 3 – End-to-End TCGA Genomics & Survival Modeling

##  Objective

Perform complete end-to-end genomic data analysis on TCGA UCEC dataset (391 patients).

---

##  Full Analytical Pipeline

###  Preprocessing
- Expression filtering (90th quantile threshold)
- Outlier detection
- Clinical data alignment
- Batch effect inspection

---

### Feature Selection
- Identified mRNA biomarkers associated with POLE mutation status
- Compared clustering before and after feature selection

---

###  Classification
- Built predictive models for mutation status
- Evaluated with ROC and confusion matrices
- Compared classifier families

---

###  Statistical Testing
- Tested association between age and selected mRNAs
- Evaluated stage vs mutation status
- Verified assumptions and used alternative tests where necessary

---

### Survival Analysis
- Kaplan–Meier survival curves
- Censoring variable engineering
- MSI vs POLE mutation outcome comparison

---

##  Skills Demonstrated

- High-dimensional feature reduction
- Model validation
- Clinical data integration
- Survival modeling
- Statistical reasoning and interpretation

---

## Why This Project Stands Out

This project represents a full machine learning workflow applied to real TCGA cancer genomics data — from preprocessing to predictive modeling to clinical outcome interpretation.