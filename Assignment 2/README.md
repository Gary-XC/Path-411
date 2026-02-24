# Mini Assignment 2 – RNA-Seq Alignment & Differential Expression Pipeline

## Objective

Build a complete RNA-seq alignment and quantification pipeline using the HISAT2–StringTie–Ballgown workflow on a high-performance computing cluster.

---

## Dataset

RNA-seq FASTQ files for:
- UHR cells
- HBR cells

Reference genome:
- Human GRCh38 assembly

---

## End-to-End Pipeline

### Quality Control
- Generated FastQC reports
- Evaluated base quality, GC content, duplication rates
- Assessed need for trimming or filtering

---

### Genome Indexing
- Built HISAT2 index from GRCh38 reference genome
- Managed memory allocation on cluster environment

---

### Read Alignment
- Aligned paired-end reads using HISAT2
- Generated SAM files
- Converted and sorted to BAM using samtools

---

### Transcript Assembly & Quantification
- Assembled transcripts using StringTie
- Generated:
  - Raw count matrices
  - TPM and FPKM normalized counts
- Merged transcript models across samples

---

### Differential Expression Preparation
- Combined count matrices in MATLAB
- Generated gene-level and transcript-level count matrices
- Prepared dataset for downstream statistical modeling

---

## Skills Demonstrated

- High-performance computing (HPC)
- Shell scripting & reproducibility
- RNA-seq alignment fundamentals
- Transcript quantification
- Multi-tool pipeline integration
- Memory and file management at scale

---

## ML Relevance

This project demonstrates:

- Large-scale raw data processing
- Pipeline engineering
- Multi-stage transformation workflows
- Reproducible computational systems

These are core competencies in production ML environments.