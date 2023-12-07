##########################
# This is a helper function for running decontX in celda package


decX <- function(seu){
  
  
  ##########################
  # This was written in Oct 6th, 2022 by MJB
  # USAGE: decontX(seu)
  #
  # seu <- pre-filtered single sample Seurat object
  ##########################
  
  require(Seurat)
  require(celda)
  
  # Convert Seurat obj to SCE, run decontX, and return Seurat obj
  counts <- GetAssayData(object = seu, slot = "counts")
  sce <- SingleCellExperiment(list(counts = counts))
  sce <- decontX(sce)
  seu[["decontXcounts"]] <- CreateAssayObject(counts = decontXcounts(sce))
  
  return(seu)
}