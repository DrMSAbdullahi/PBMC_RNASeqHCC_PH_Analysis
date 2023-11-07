# Persistent homology identifies pathways associated with HCC from PBMC samples

This is a supplementary material for codes used in the paper:
Abdullahi et al., "Persistent homology identifies pathways associated with hepatocellular carcinoma from peripheral blood samples", Submitted to "Mathematics" November, 2023.

## Table of Contents

- [Introduction](#introduction)
- [Datasets](#datasets)
- [Topological Data Analysis](#topological-data-analysis)
- [Permutation Test](#permutation-test)
- [Pathway Enrichment Analysis](#pathway_enrichment-analysis)
- [License](#license)

## Introduction

This project is designed to complement the codebase used in our research paper. In the paper, we employed topological data analysis, specifically the persistent homology method, to identify crucial pathways associated with hepatocellular carcinoma using peripheral blood mononuclear cell (PBMC) samples.

## Datasets

We used the following datasets in our project. 

- **RNA-Seq Gene Expression Dataset**: This dataset comprises 26,575 gene expressions from 34 samples, with 17 representing healthy samples and 17 representing disease samples. The original dataset was sourced from the study by Han et al.[1], titled 'RNA-seq profiling reveals PBMC RNA as a potential biomarker for hepatocellular carcinoma.' For access to the original sequencing data, please refer to NCBI project PRJNA739257 [here](https://dataview.ncbi.nlm.nih.gov/object/PRJNA739257).

- **KEGG Pathways Dataset**: We obtained KEGG pathway datasets from the KEGG database using R software (version 4.3.0). Our focus was specifically on the 230 signaling and disease pathways, each of which consists of a minimum of 25 genes. Additionally, at least 70% of the genes from each pathway were required to be present in our gene expression dataset.

## Topological Data Analysis

### Project Dependencies

- **Gudhi Library (Python Version)**: We utilized the 'Gudhi' library [2] (version 3.7.1) for performing computations and visualizations of persistent homology. You can find the Python version of 'Gudhi' at [https://gudhi.inria.fr/](https://gudhi.inria.fr/).

### Persistent Homology Computation

- **Inter-Sample Dissimilarity Distance Matrix**: The persistent homology computations were based on the inter-sample dissimilarity distance matrix, calculated using the complement of the Pearson correlation coefficient (1 - p).

- **VR Complexes Construction**: We adopted the VR complexes construction method to create the simplicial complexes for our sample point clouds.

### Topological Descriptors

- **User-Friendly Functions**: To simplify the calculation of topological descriptors, we have provided user-friendly functions that assist with the computations and further analysis.

### Visualization

- **Persistence Diagrams and Barcodes**: All persistence diagrams and persistence barcodes were respectively plotted using the 'plot_persistence_diagram' and 'plot_persistence_barcode' functions from the 'Gudhi' library in Python.

- **Betti Curves**: We demonstrated the Betti curves using our defined functions.


```bash
# Example installation instructions
npm install
npm start
```

## Permutation Test

Two-tailed statistical permutation tests was performed using Python (version 3.10.7) and the Benjamini–Hochberg correction method was adopted for multiple testing using 'statsmodels' package (version 0.13.5). The \ac{fdr} threshold was set to less than 0.05 throughout.

## Gene Diffrential Expression Analysis

Differential expression analysis was performed using the `PyDESeq2' \cite{muzellec2023pydeseq2} package (version 0.3.0) in Python (version 3.10.7). Only genes with adjusted p-value $< 0.05$ (adjusted using the Benjamini-Hochberg method), and log fold change $\geq 1$ were considered as significantly differentially expressed.

## Pathway Enrichment Analysis

Differential enrichment analysis of pathways was carried out using `Scipy' package (version 1.10.1) in Python (version 3.10.7).

## License

Statement on License goes here.



[1] Han, Z. et al. (2021). 'RNA-seq profiling reveals PBMC RNA as a potential biomarker for hepatocellular carcinoma.' Sci. Reports, 11, 17797.

[2] Maria, C., Oudot, S.Y., & Glisse, M. (2014). 'Gudhi: Merging computational topology and geometric data analysis.' Proceedings of the 30th International Symposium on Computational Geometry (SoCG).
