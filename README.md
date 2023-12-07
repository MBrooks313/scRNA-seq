# scRNA-seq
Single cell RNA-seq analysis pipeline.

This repo has a collection of scripts to process scRNA-seq data. There is a kallisto-bustools version and a CellRanger pipeline for primary data processing.

**Kalliso-bustools** 
[kallisto-bustools](https://www.kallistobus.tools/) output to  

**CellRanger**
[CellRanger](https://www.10xgenomics.com/support/software/cell-ranger/getting-started/cr-what-is-cell-ranger)


The ##src## directory contains scripts for kallisto or CellRanger import in Seurat object. 



Both primary pipelines can be imported in R for creation of a SeuraT OBJECT. [Seurat](https://satijalab.org/seurat/articles/get_started.html)

[DecontX](https://bioconductor.org/packages/release/bioc/vignettes/celda/inst/doc/decontX.html#running-decontx) removes ambiant RNA contamination.

[DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder) performs automated detection of heterotypic doublets. 

[SCT v2](https://satijalab.org/seurat/articles/sctransform_v2_vignette.html) is used for normalization.

[scType](https://github.com/IanevskiAleksandr/sc-type) is used for cell type annotation with genes determined with Sanes' laboratory scRNA-seq data. 


