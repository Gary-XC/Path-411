# Mini Assignment 1 – Genomic Data Preprocessing Utilities

## Objective

Developing reusable MATLAB functions to support preprocessing and quality control of high-dimensional genomic datasets.

---

## Implemented Functions

### Relative Frequency Normalization
- Implemented column-wise normalization for genomic count matrices
- Converted raw counts into relative frequencies
- Ensured numerical stability and avoided inefficient looping

---

### Sample Correlation Analysis
- Computed average Spearman or Pearson correlation per sample
- Implemented vectorized solution (no unnecessary loops)
- Enabled detection of outlier samples via correlation structure

---

### Low Expression Feature Filtering
- Implemented quantile-based filtering function
- Identified low-expression miRNA features
- Reduced dimensionality prior to modeling

---

## Skills Demonstrated

- Vectorized numerical computing
- Matrix manipulation at scale
- Reusable function design
- Defensive programming practices
- Genomic data normalization principles

---

## Relevance to ML

This project mirrors real-world ML preprocessing pipelines:

- Feature filtering
- Noise reduction
- Outlier detection
- Data conditioning prior to modeling