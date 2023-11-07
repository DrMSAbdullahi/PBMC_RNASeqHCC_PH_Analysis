# Persistent homology identifies pathways associated with HCC from PBMC samples

This is a supplementary material for codes used in the paper:
Abdullahi et al., "Persistent homology identifies pathways associated with hepatocellular carcinoma from peripheral blood samples", Submitted to "Mathematics" November, 2023.

## Table of Contents

- [Introduction](#introduction)
- [Datasets](#datasets)
- [Topological Data Analysis](#topological-data-analysis)
- [Permutation Test](#permutation test)
- [Pathway Enrichment Analysis](#pathway_enrichment-analysis)
- [License](#license)

## Introduction

This project is designed to complement the codebase used in our research paper. In the paper, we employed topological data analysis, specifically the persistent homology method, to identify crucial pathways associated with hepatocellular carcinoma using peripheral blood mononuclear cell (PBMC) samples.

## Datasets

We used the following datasets in our project. 

- **RNA-Seq Gene Expression Dataset**: This dataset comprises 262,575 gene expressions from 34 samples, with 17 representing healthy samples and 17 representing disease samples. The original dataset was sourced from the study by Han et al., titled 'RNA-seq profiling reveals PBMC RNA as a potential biomarker for hepatocellular carcinoma.' For access to the original sequencing data, please refer to NCBI project PRJNA739257 [here](https://dataview.ncbi.nlm.nih.gov/object/PRJNA739257).

- **KEGG Pathways Dataset**: We obtained KEGG pathway datasets from the KEGG database using R software (version 4.3.0). Our focus was specifically on the 230 signaling and disease pathways, each of which consists of a minimum of 25 genes. Additionally, at least 70% of the genes from each pathway were required to be present in our gene expression dataset.

## Topological Data Analysis

Include instructions on how to get started with your project. This should include installation steps, prerequisites, and any initial setup.

```bash
# Example installation instructions
npm install
npm start
```

## Permutation Test

Provide a brief and clear introduction to your project. Explain what it does and why it's useful.

## Pathway Enrichment Analysis

Provide a brief and clear introduction to your project. Explain what it does and why it's useful.

## License

Statement on License goes here.
