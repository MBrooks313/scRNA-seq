##########################
# This function imports CellRanger scRNA-seq h5 files into a Seurat object, runs decontX and doublet finder and exports Seurat list of samples.



cr_import <- function(files, expt_nm, min.cells=5, min.features=200, max.mt=20, dX=TRUE, dF=TRUE, dub_exp=0.075, dub_counts="RNA", dub_remove=FALSE, helper_dir="../src/data", out_dir="../data/processed"){
  
  
  
  ##########################
  # This was written in Dec 1st, 2022 by MJB
  # USAGE: cr_import(run_info_files, expt_nm, min.cell = 5, min.features = 200, max.mt = 20, helper_dir = "../src/data", out_dir = "../data/processed")
  #
  # run_info_files <- named vector of paths to the sample kallisto run_info.json
  # expt_nm <- name of the experiment, will be name of output directory of files
  # min.cells <- minumum number of cells for gene to be kept, default = 5 
  # min.features <- minimum genes per cell to be kept, default = 200
  # max.mt <- maximum mitochondrial content per cell
  # dX <- run the decontX protocol
  # dF <- run the DoubletFinder protocol
  # dub_exp <- expected percentage of doublets (based on 10X protocol manual)
  # dub_counts <- count assay to use in doubletFinder, one of c("RNA", "decontXcounts")
  # dub_remove <- whether to remove the cells maker as doublets (default = FALSE)
  # helper_dir <- directory of helper scripts (tr2g.R, read_count_output.R)
  # out_dir <- directory of output directories (tr2g.R, read_count_output.R)
  #
  # Value
  # seu_list - Seurat object list with doublet cells marked/filtered out
  # dub_df - data.frame of resulting cell counts at each filtration step
  ##########################
  
  
  require(tidyverse)
  require(Seurat)
  require(DoubletFinder)
  
  source(file.path(helper_dir, "doublet_Finder.R"))
  source(file.path(helper_dir, "decontX.R"))
  
  # Empty lists for return
  seu_list <- list()
  dub_df <- data.frame()
  
  # Prepare output directories
  out_expt <- file.path(out_dir, expt_nm)
  out_qc <- file.path(out_dir, expt_nm, "QC_indi")
  
  dir.create(out_qc, recursive = T)
  
  
  # Loop for each sample in run_info_list
  for (i in 1:length(files)){
    
    samp <- files[i]
    samp_nm <- names(files)[i]
    print(paste("Importing matrix file: ", samp_nm))

    # Import 10x h5 object
    cr_h5 <- Read10X_h5(samp, use.names = TRUE, unique.features = TRUE)
    
    # Create Seurat object
    print("Creating Seurat object...")
    seu_tmp <- CreateSeuratObject(cr_h5, min.cells = min.cells, min.features = min.features, project = samp_nm)
    Idents(seu_tmp) <- samp_nm
    
    filt_min <- dim(seu_tmp)[2]
    
    # Calculate the miochondria content per cell
    print("Filtering for mito content...")
    seu_tmp[["percent.mt"]] <- PercentageFeatureSet(seu_tmp, pattern = "^mt-")
    
    # Subset for cells with less than 20% MT
    seu_tmp <- subset(seu_tmp, subset = percent.mt < max.mt)
    mt_cells <- dim(seu_tmp)[2]


    # Run decontX to get rid of ambiant RNA contamination
    if (dX){
      
      print("Running decontX...")
      seu_tmp <- decX(seu_tmp)
    }
    

    # Run DoubletFinder
    if (dF){
      
      print("Running DoubletFinder...")
    
      seu_tmp <- dub_find(seu_tmp, expect_dub = dub_exp, sct = TRUE, assay_counts = dub_counts)
    
      # Plots for doublet finder
      print("Making DoubletFinder plots...")

      df_nm <- tail(colnames(seu_tmp@meta.data), n = 1)
      a <- DimPlot(seu_tmp, group.by = df_nm)
      b <- VlnPlot(seu_tmp, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, log = T, group.by = df_nm)
      
      cowplot::ggsave2(plot = a, file = file.path(out_qc, paste0(expt_nm, "_Dim_Dub_", samp_nm, ".png")), width = 12, height = 8)
      cowplot::ggsave2(plot = b, file = file.path(out_qc, paste0(expt_nm, "_Viol_feat-count_", samp_nm, ".png")), width = 6, height = 4)

      # Make a Seurat object for filtered cells
      seu_tmp <- AddMetaData(seu_tmp, seu_tmp@meta.data[[df_nm]], col.name = "DubFind")
      
      # Filter out doublets
      if (dub_remove){
        print("Filter out doublets...")
        seu_tmp <- subset(seu_tmp, subset = DubFind == "Singlet")
      }
      
    }
    
    
    # Get cell numbers for data.frame
    og_cells <- dim(cr_h5)[2]
    filt_cells <- sum(seu_tmp$DubFind == "Singlet")
    
    tmp_df <- tibble(c("orig_cells", 
                       "orig_filt", 
                       "mt_filt", 
                       "dub_filt",
                       "filt_cells"),
                     c(og_cells, 
                       og_cells - filt_min, 
                       filt_min - mt_cells,
                       mt_cells - filt_cells,
                       filt_cells))
    tmp_df$samp <- samp_nm
    dub_df <- rbind(dub_df, tmp_df)
    
    # Add Seurat object to list
    seu_list[[samp_nm]] <- seu_tmp
    
  }
  
  # Return objects
  seu_out <- list("dub_df" = dub_df,
                  "seu_list" = seu_list)
        
  return(seu_out)
  
}


