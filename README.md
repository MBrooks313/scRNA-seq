## scRNA-seq
Single cell RNA-seq analysis pipeline Kallisto/Bustools-Seurat.

This repo contains scripts for processing [kallisto-bustools](https://www.kallistobus.tools/) output to [Seurat](https://satijalab.org/seurat/articles/get_started.html) object. [DecontX](https://bioconductor.org/packages/release/bioc/vignettes/celda/inst/doc/decontX.html#running-decontx) and [DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder) are used to remove ambiant RNA contamination and automated detection of heterotypic doublets. [SCT v2](https://satijalab.org/seurat/articles/sctransform_v2_vignette.html) is used for normalization and [scType](https://github.com/IanevskiAleksandr/sc-type) is used for cell type annotation with genes determined with Sanes' laboratory scRNA-seq data. 

### kallisto/bustools

Kallisto/bustools scripts are run on an HPC using Slurm Workload Manager.

00_config.sh loads the python 3 environment containing the kallisto/bustools software.  