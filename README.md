# KODAMA
An unsupervised and semi-supervised learning algorithm to performs feature extraction from noisy and high-dimensional data

## News 

KODAMA facilitates identification of patterns representing underlying groups on all samples in a data set. 
This is an improved version of KODAMA algorithm for spatially-aware dimensionality reduction. A landmarks procudere has been implemented to adapt the algorithm to the analysis of data set with more than 10,000 entries. 

The KODAMA package has been integrated with t-SNE and UMAP to convert the KODAMA's dissimilarity matrix in a low dimensional space. 



[Zinga, M. M., Abdel-Shafy, E., Melak, T., Vignoli, A., Piazza, S., Zerbini, L. F., ... & Cacciatore, S. (2022). KODAMA exploratory analysis in metabolic phenotyping. Frontiers in Molecular Biosciences, 9.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9887019/)

[Cacciatore, S., Tenori, L., Luchinat, C., Bennett, P. R., & MacIntyre, D. A. (2017). KODAMA: an R package for knowledge discovery and data mining. Bioinformatics, 33(4), 621-623.](https://academic.oup.com/bioinformatics/article/33/4/621/2667156?login=false)

[Cacciatore, S., Luchinat, C., & Tenori, L. (2014). Knowledge discovery by accuracy maximization. Proceedings of the National Academy of Sciences, 111(14), 5117-5122.](https://www.pnas.org/doi/abs/10.1073/pnas.1220873111)



## Installation

The KODAMA is avialable on https://CRAN.R-project.org/package=KODAMA.

```
library(devtools)
install_github("tkcaccia/KODAMA")

```


## Applications 
Here below, we introduced three different applications of the KODAMA algorithm.

1.  [Metabolomic data](https://github.com/tkcaccia/KODAMA/blob/main/docs/Metabolomics_data.md).

2.  [Single cell RNA seq data](https://github.com/tkcaccia/KODAMA/blob/main/docs/Single_cell_RNA_seq.md).

3.  [Spatial Transcriptomic data](https://github.com/tkcaccia/KODAMA/blob/main/docs/Spatial%20_transcriptomic.md).
