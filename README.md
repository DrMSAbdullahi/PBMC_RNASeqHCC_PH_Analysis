# Persistent homology identifies pathways associated with HCC from PBMC samples

This is a supplementary material for codes used in the paper:
Abdullahi et al., "Persistent homology identifies pathways associated with hepatocellular carcinoma from peripheral blood samples", Mathematics 2024, 12(5), 725.

## Table of Contents

- [Introduction](#introduction)
- [Datasets](#datasets)
- [Topological Data Analysis](#topological-data-analysis)
- [Permutation Test](#permutation-test)
- [Pathway Enrichment Analysis](#pathway_enrichment-analysis)
- [License](#license)

## Introduction

This project is designed to complement the codebase used in our research paper. In the paper, we employed topological data analysis, specifically the persistent homology method, to identify key pathways for early detection of hepatocellular carcinoma (HCC) from peripheral blood mononuclear cell (PBMC) samples.

## Datasets

We used the following datasets in our project. 

- **RNA-Seq Gene Expression Dataset**: This dataset comprises 26,575 gene expressions from 34 samples, with 17 representing healthy samples and 17 representing disease samples. The original dataset was sourced from the study by Han et al.[1], titled 'RNA-seq profiling reveals PBMC RNA as a potential biomarker for hepatocellular carcinoma.' To access the original sequencing data, please refer to NCBI project PRJNA739257 [here](https://dataview.ncbi.nlm.nih.gov/object/PRJNA739257).

- **KEGG Pathways Dataset**: We obtained KEGG pathway datasets from the KEGG database using R software (version 4.3.0). Our focus was specifically on the 251 metabolic and signaling pathways.

## Topological Data Analysis

### Project Dependencies

- **Gudhi Library (Python Version)**: We utilized the 'Gudhi' library [2] (version 3.7.1) for performing computations and visualizations of persistent homology. You can find the Python version of 'Gudhi' at [https://gudhi.inria.fr/](https://gudhi.inria.fr/).

### Persistent Homology Computation

- **Inter-Sample Dissimilarity Distance Matrix**: The persistent homology computations were based on the inter-sample dissimilarity distance matrix, calculated using the complement of the Pearson correlation coefficient p (i.e., 1 - p).

- **VR Complexes Construction**: We adopted the VR complexes construction method to create the simplicial complexes for our sample point clouds.

### Topological Descriptors

- **User-Friendly Functions**: To simplify the calculation of topological descriptors, we have provided user-friendly functions that assist with the computations and further analysis.

### Visualization

- **Persistence Diagrams and Barcodes**: All persistence diagrams and persistence barcodes were plotted using the 'plot_persistence_diagram' and 'plot_persistence_barcode' functions from the 'Gudhi' library in Python.

- **Betti Curves**: We demonstrated the Betti curves using our own defined functions.

- **Permutation Test Plots**: We employed the `seaborn.kdeplot` function from the Seaborn data visualization package for generating KDE (Kernel Density Estimation) plots of our permutation distributions. This allowed us to visualize the density estimators for our permutation distributions.

## Permutation Test

A two-tailed statistical permutation test was performed using Python (version 3.10.7) and the Benjamini–Hochberg correction method was adopted for multiple testing using 'statsmodels' package (version 0.13.5). The False Discovery Rate (FDR) threshold was set to less than 0.05 throughout.

## Gene Differential Expression Analysis

Differential expression analysis was performed using the `PyDESeq2' [3] package (version 0.4.4) in Python (version 3.10.7). Only genes with adjusted p-value $< 0.05$ (adjusted using the Benjamini-Hochberg method), and log fold change $\geq 1$ were considered as significantly differentially expressed.

## Pathway Enrichment Analysis

Differential enrichment analysis of pathways was carried out using the `Scipy' package (version 1.10.1) in Python (version 3.10.7).

## License

This codebook is licensed under the [Creative Commons Attribution 4.0 International License](https://creativecommons.org/licenses/by/4.0/).

**You are free to:**
- **Share:** Copy and redistribute the material in any medium or format.
- **Adapt:** Remix, transform, and build upon the material for any purpose, even commercially.

**Under the following terms:**
- **Attribution:** You must give appropriate credit, provide a link to the license, and indicate if changes were made. You may do so in any reasonable manner, but not in any way that suggests the licensor endorses you or your use.

**No additional restrictions:** You may not apply legal terms or technological measures that legally restrict others from doing anything the license permits.

This codebook is provided "as is," without any warranty or guarantee of fitness for a particular purpose. The codebook's users are solely responsible for their use and any consequences thereof.

If you use this codebook in your research or publication, please cite the original paper:

[Abdullahi, M.S.; Suratanee, A.; Piro, R.M.; Plaimas, K. "Persistent Homology Identifies Pathways Associated with Hepatocellular Carcinoma from Peripheral Blood Samples". Mathematics 2024, 12, 725. https://doi.org/10.3390/math12050725]

For questions or additional permissions, contact [abdullahi.sirajo@udusok.edu.ng].

#### References

[1] Han, Z., Feng, W., Hu, R., Ge, Q., Ma, W., Zhang, W., Xu, S., Zhan, B., Zhang, L., Sun, X. et al. (2021). 'RNA-seq profiling reveals PBMC RNA as a potential biomarker for hepatocellular carcinoma.' Sci. Reports, 11, 17797.

[2] Maria, C., Oudot, S.Y., Glisse, M. (2014). 'Gudhi: Merging computational topology and geometric data analysis.' Proceedings of the 30th International Symposium on Computational Geometry (SoCG).

[3] Muzellec, B., Teleńczuk, M., Cabeli, V., Andreux, M. (2023). 'PyDESeq2: A Python Package for Bulk RNA-Seq Differential Expression Analysis.' Bioinformatics, 39, btad547.
