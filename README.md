# Response Specificity
Quantifying the Response Specificity of Macrophages (Sheu et al, 2023)

## In Brief
Emerging evidence suggests that macrophages are capable of immune threat-appropriate responses. To quantify response specificity, we introduce experimental and analytical workflows. Our studies reveal alterations by microenvironmental cytokines in vitro and within unhealthy mice, suggesting a means to assess the functional state of macrophages.

## Outline
This repository contains code for generating computational models to quantify macrophage Response Specificity. 
The input counts matrix containing single cell RNA abundance log normalized counts (see Methods) can be found in the output subfolder, 
along with correspondence of cell barcodes to meta data.
The script figures_run.R contains compiled code from all figures in the paper. In particular,
- code for training machine learning models that can be used to assess ligand specificity (Figure 1).
- code for information theoretic frameworks that are employed to quantify the Response Specificity of genes and genesets (Figure 4).
- code for identifying genes that contribute to a summary metric for Response Specificity, called the Response Specificity Index (Figure 6). 

Intermediate files that are generated as the result of more time-intensive computations are in the folders ./infotheo and ./ml_models.
Processed scRNAseq data related to Figure 1 (naive macrophages), Figure 5 (polarized macrophages), and Figure 7 (peritoneal macrophages), are provided as .rds files in the folder ./output. 

Machine learning and information theory calculations were run on Dual CPU Intel Xeon E5-2690 v3 @ 2.60GHz. 

## Citation
Sheu, Katherine M., Aditya A. Guru, and Alexander Hoffmann. “Quantifying Stimulus-Response Specificity to Probe the Functional State of Macrophages.” Cell Systems 14, no. 3 (March 15, 2023): 180-195.e5. https://doi.org/10.1016/j.cels.2022.12.012.
