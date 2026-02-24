# Group Assignment 2 – Feature Selection & Classification

##  Objective

Identify discriminative genomic features and evaluate classification performance on GI-NET cancer dataset.

---

## Dataset

Normalized miRNA expression data  
4 tumor types:
- Pancreatic
- Ileal
- Appendiceal
- Rectal

---

## Unsupervised Learning

- t-SNE visualization of class separability
- Hierarchical clustering (multiple linkage methods)
- Median centering and log transformations

---

##  Feature Selection Methods

Explored multiple approaches:
- Chi-square
- mRMR
- ReliefF
- Random forest importance
- Laplacian score

Compared importance scores and selected optimal subsets.

---

##  Supervised Learning

Benchmarked classifiers:
- Logistic regression
- SVM
- KNN
- Decision trees
- Ensemble models

Evaluated with:
- Cross-validation
- ROC curves
- Confusion matrices

---

##  Key Insight

Demonstrated how feature selection significantly improves separability and model performance in high-dimensional data.