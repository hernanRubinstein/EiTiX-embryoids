# Mouse-embryo model derived exclusively from embryonic stem cells undergoes neurulation and heart development

This repository is the accompanying code for our paper on the single-cell single structure analysis of Day 6 and 8 EiTiX embryoids.
The code is splitted into one jupyter notebook (to generate all figures in this work, named 'Analysis') in the main directory. 
Additional functions and data construction R scripts can be found in the 'scripts' folder (where the order of execution is assigned by a capital letter A-H).

## Raw and processed data storage
Prior to the analysis, after cloning the repository, please download first the necessary data from the following link:

[Main UMI matrices](https://eitix-embryoids.s3.eu-west-1.amazonaws.com/umi_matrices.tar.gz)

[Additional data](https://eitix-embryoids.s3.eu-west-1.amazonaws.com/scrna_db.tar.gz)

*This will download all the necessary data (UMI matrices of the scRNA-seq data used in the paper) including processed metacell objects necessary for generating the figures.*

## Necessary R packages
The analysis was done using R 4.0.5 and the following packages:

- devtools_2.4.2
- usethis_2.0.1
- here_1.0.1         
- slanter_0.2-0
- DoubletFinder_2.0.3 
- SeuratObject_4.0.2 
- Seurat_4.0.3
- forcats_0.5.1
- stringr_1.4.0      
- dplyr_1.0.9
- purrr_0.3.4
- readr_2.1.0        
- tidyr_1.2.0
- tibble_3.1.3
- ggplot2_3.3.5      
- tidyverse_1.3.1
- umap_0.2.7.0
- tgutil_0.1.13      
- tgstat_2.3.17
- metacell_0.3.7
- Matrix_1.3-4 
- data.table_1.14.2
- qvalue_2.22.0
- princurve_2.1.6
- RColorBrewer_1.1-2
- tglkmeans_0.3.4
- zoo_1.8-9
- ggrepel_0.9.1
