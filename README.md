# scRNA-seq
Single cell RNA-seq analysis pipeline.

This repo has a collection of scripts to process scRNA-seq data. There is a kallisto-bustools version and a CellRanger pipeline for primary data processing.


[**kallisto-bustools**](https://www.kallistobus.tools/) processes reads only to the annotated exons of the genes. 

[**CellRanger**](https://www.10xgenomics.com/support/software/cell-ranger/getting-started/cr-what-is-cell-ranger) processes reads to the introns as well as the exons by default. Exon only analysis can be specified in the options.

Scripts for import of both primary analysis pipelines into R for creation of a [Seurat](https://satijalab.org/seurat/articles/get_started.html) object are located in the pipelines respective directory.  

The **src** directory contains scripts for further processing the Seurat object. 

[DecontX](https://bioconductor.org/packages/release/bioc/vignettes/celda/inst/doc/decontX.html#running-decontx) removes ambiant RNA contamination.

[DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder) performs automated detection of heterotypic doublets. 

[SCT v2](https://satijalab.org/seurat/articles/sctransform_v2_vignette.html) is used for normalization.

[scType](https://github.com/IanevskiAleksandr/sc-type) is used for cell type annotation.


The **annotation** directory contains cell type annotation lists for running with scType.

 




