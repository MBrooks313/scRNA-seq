##########################
# This is a helper function for running DoubletFinder


dub_find <- function(seu, expect_dub = 0.075, sct = TRUE, assay_counts = "RNA"){
  
  
  ##########################
  # This was written in Oct 5th, 2022 by MJB
  # USAGE: dub_find(seu, expect_dub = 0.075, sct = TRUE)
  #
  # seu <- pre-filtered single sample Seurat object
  # expect_dub <- expected doublet rate based on single cell density used in experiment
  # sct <- was SCTransform used for normalization
  ##########################
  
  
  require(Seurat)
  require(tidyverse)
  require(DoubletFinder)
  
  # Prep Seurat object
  seu <- SCTransform(seu, assay = assay_counts) %>% 
    RunPCA() %>%
    RunUMAP(dims = 1:10)
  
  ## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
  print("Doublet finder...")
  sweep.res.list <- paramSweep_v3(seu, PCs = 1:10, sct = sct)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK_val <- as.numeric(as.character(bcmvn$pK[which(bcmvn$BCmetric == max(bcmvn$BCmetric))]))
  
  ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
  homotypic.prop <- modelHomotypic(seu@meta.data$seurat_clusters)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  nExp_poi <- round(expect_dub*nrow(seu@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
  seu_dub <- doubletFinder_v3(seu, PCs = 1:10, pN = 0.25, pK = pK_val, nExp = nExp_poi, reuse.pANN = FALSE, sct = sct)
  pANN_nm <- paste0("pANN_0.25_", pK_val, "_", nExp_poi)
  seu_dub <- doubletFinder_v3(seu_dub, PCs = 1:10, pN = 0.25, pK = pK_val, nExp = nExp_poi.adj, reuse.pANN = pANN_nm, sct = sct)
  
  ## Return
  return(seu_dub)
  
}



