# ResponseSpecificity
Quantifying the Response Specificity of macrophages

This repository contains code for generating computational models to quantify macrophage Response Specificity. 
The input counts matrix containing single cell RNA abundance log normalized counts (see Methods) can be found in the output subfolder, 
along with correspondence of cell barcodes to meta data.
The script figures_run.R contains compiled code from all figures in the paper. In particular,
- code for training machine learning models that can be used to assess ligand specificity (Figure 1).
- code for information theoretic frameworks that are employed to quantify the Response Specificity of genes and genesets (Figure 4).
- code for identifying genes that contribute to a summary metric for Response Specificity, called the Response Specificity Score (Figure 6). 

Intermediate files that are generated as the result of more time-intensive computations are in the folders ./infotheo and ./ml_models.
Processed scRNAseq data related to Figure 1 (naive macrophages), Figure 5 (polarized macrophages), and Figure 7 (peritoneal macrophages),are provided as .rds files in the folder ./output. 

Machine learning and information theory calculations were run on Dual CPU Intel Xeon E5-2690 v3 @ 2.60GHz. 