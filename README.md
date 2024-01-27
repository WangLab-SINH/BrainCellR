# BrainCellR: A Precise Cell Type Nomenclature R Package for Comparative Analysis Across Brain Single-Cell Datasets  

`BrainCellR` provides researchers with a powerful and user-friendly package for efficient cell type classfication and nomination of single-cell transcriptomic data in the brain.

`BrainCellR` process comprises the following steps:  

* Preparation of input data
* Fine-scale clustering of single-cell data
* Provision of clusters to a supervised cell type classifier capable of outputting major cell classes and subclasses
* Identification of differentially expressed genes and select marker genes from the identified differentially expressed genes
  + Find differentially expressed genes  
  + Selection of genes that exhibit a high fold change value and are expressed in the majority of cells within a given cell type
  + Ranking of genes based on specific score   
* Acquisition of the final cell type annotation
* Identification of markers that cannot be detected using the ROC method

  
## Installation

`BrainCellR` has a dependencie from Github:
```
devtools::install_github("AllenInstitute/scrattch.hicat")
```
Also, WGCNA and SingleR

`BrainCellR` can be installed with:
```
devtools::install_github("WangLab-SINH/BrainCellR")
```
## Tutorials

Tutorials can be seen at https://wanglab-sinh.github.io/braincellr<br>
Data used in tutorial can be downloaded from https://github.com/WangLab-SINH/WangLab-SINH.github.io

Reference data for cell type class and subclass classification can be obtained from https://drive.google.com/drive/folders/1q9JT0JFhBvc6CvkbzZVXmEANXfgQ8VID?usp=sharing

You can contact chiyuhao2018@sinh.ac.cn for any questions.
