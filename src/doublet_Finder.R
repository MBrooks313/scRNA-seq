##########################
# This is a helper function for running DoubletFinder


dub_find <- function(seu, expect_dub = 0.075, sct = TRUE, assay_counts = "RNA", vst = "v2", PC_num = NULL,  helper_dir="../src/data"){
  
  
  ##########################
  # This was written in Oct 5th, 2022 by MJB
  # USAGE: dub_find(seu, expect_dub = 0.075, sct = TRUE)
  #
  # seu <- pre-filtered single sample Seurat object
  # expect_dub <- expected doublet rate based on single cell density used in experiment
  # sct <- was SCTransform used for normalization
  # vst <- version of SCT, NULL = v1, v2 = v2
  # PC_num <- specify number of PCs to use in DoubletFinder, NULL = optimal PC
  # helper_dir <- location of helper files
  ##########################
  
  
  require(Seurat)
  require(tidyverse)
  require(DoubletFinder)
  source(file.path(helper_dir, "opt_PC.R"))
  
  
  # Find optimal PC
  if (is.null(PC_num)){
    seu <- SCTransform(seu, assay = assay_counts, vst.flavor=vst) %>% 
      RunPCA() %>%
      RunUMAP(dims = 1:30)
    pc <- opt_PC(seu)
  } else{
    pc <- PC_num
  }
  
  
  # Prep Seurat object
  seu <- SCTransform(seu, assay = assay_counts, vst.flavor=vst) %>% 
    RunPCA() %>%
    RunUMAP(dims = 1:pc)
  
  
  ## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
  print("Doublet finder...")
  sweep.res.list <- paramSweep_v3(seu, PCs = 1:pc, sct = sct)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK_val <- as.numeric(as.character(bcmvn$pK[which(bcmvn$BCmetric == max(bcmvn$BCmetric))]))
  
  ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
  homotypic.prop <- modelHomotypic(seu@meta.data$seurat_clusters)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  nExp_poi <- round(expect_dub*nrow(seu@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
  seu_dub <- doubletFinder_v3(seu, PCs = 1:pc, pN = 0.25, pK = pK_val, nExp = nExp_poi, reuse.pANN = FALSE, sct = sct)
  pANN_nm <- paste0("pANN_0.25_", pK_val, "_", nExp_poi)
  seu_dub <- doubletFinder_v3(seu_dub, PCs = 1:pc, pN = 0.25, pK = pK_val, nExp = nExp_poi.adj, reuse.pANN = pANN_nm, sct = sct)
  
  ## Return
  return(seu_dub)
  
}



