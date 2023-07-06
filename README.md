# BrainCellR: A Precise Cell Type Nomenclature R Package for Comparative Analysis Across Brain Single-Cell Datasets  

`BrainCellR` provides researchers with a powerful and user-friendly package for efficient cell type nomination of single-cell transcriptomic data.

`BrainCellR` consists of the following steps:  

* Input data preparation
* Cluster single-cell data at a fine scale  
* Supply the clusters to a supervised cell type classifier which can output cell major classes and subclasses
* Find differentially expressed genes and select marker genes from the identified differentially expressed genes
  + Find differentially expressed genes  
  + Select genes that have high fold change value and expressed in most of the cells belonging to the cell type
  + Rank the genes by binary score    
* Get the final cell type annotation
* Find markers that can not be identified by ROC method

  
## Installation

`BrainCellR` has a dependencie from Github:
```
devtools::install_github("AllenInstitute/scrattch.hicat")
```

`BrainCellR` can be installed with:
```
devtools::install_github("WangLab-SINH/BrainCellR")
```
## Tutorials

Tutorials can be seen at https://wanglab-sinh.github.io/braincellr<br>
