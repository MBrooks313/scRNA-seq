##########################
# This function imports kallisto scRNA-seq files into a matrix, runs doublet finder and exports Seurat list of samples



kallisto_import <- function(files, expt_nm, bus_filt=TRUE, min.cells=5, min.features=200, max.mt=20, dX=TRUE, dF=TRUE, dub_exp=0.075, dub_counts="decontXcounts", helper_dir="../src/data", out_dir="../data/processed", gtf){
  
  
  ##########################
  # This was written in Sept 29th, 2022 by MJB
  # USAGE: kallisto_import(run_info_files, expt_nm, bus_filt = TRUE, min.cell = 5, min.features = 200, max.mt = 20, helper_dir = "../src/data", out_dir = "../data/processed")
  #
  # run_info_files <- named vector of paths to the sample kallisto run_info.json
  # expt_nm <- name of the experiment, will be name of output directory of files
  # bus_filt <- use bus filtered count, default = TRUE
  # min.cells <- minumum number of cells for gene to be kept, default = 5 
  # min.features <- minimum genes per cell to be kept, default = 200
  # max.mt <- maximum mitochondrial content per cell
  # dX <- run the decontX protocol
  # dF <- run the DoubletFinder protocol
  # dub_exp <- expected percentage of doublets (based on 10X protocol manual)
  # dub_counts <- count assay to use in doubletFinder, one of c("RNA", "decontXcounts")
  # helper_dir <- directory of helper scripts (tr2g.R, read_count_output.R)
  # out_dir <- directory of output directories (tr2g.R, read_count_output.R)
  # gtf <- location of the gtf file for transcript to gene table creation
  #
  # Value
  # seu_list - Seurat object list with doublet cells filtered out
  # dub_df - data.frame of resulting cell counts at each filtration step
  ##########################
  
  
  require(DropletUtils)
  require(tidyverse)
  require(Seurat)
  require(DoubletFinder)
  
  source(file.path(helper_dir, "get_tr2g.R"))
  source(file.path(helper_dir, "read_count_output.R"))
  source(file.path(helper_dir, "doublet_Finder.R"))
  source(file.path(helper_dir, "decontX.R"))
  
  # Empty lists for return
  seu_list <- list()
  dub_df <- data.frame()
  
  # Prepare output directories
  out_expt <- file.path(out_dir, expt_nm)
  out_qc <- file.path(out_dir, expt_nm, "QC_indi")
  
  dir.create(out_qc, recursive = T)
  
  # Make tr2g data.frame
  print("Getting transcript to gene information...")
  tr2g <- get_tr2g(gtf)
  
  # Loop for each sample in run_info_list
  for (i in 1:length(files)){
    
    samp <- files[i]
    samp_nm <- names(files)[i]
    print(paste("Importing matrix file: ", samp_nm))
    
    # Import data using bustools filtered data or not
    ifelse(bus_filt, bus <- "counts_filtered", bus <- "counts_unfiltered")
    
    dir_filt <- gsub("/run_info.json", "", samp)
    dir_filt <- paste(dir_filt, bus, sep = "/")
    res_mat_filt <- read_count_output(dir_filt[1])
    
    # Prep Seurat object
    # Subset matrix for only gene IDs in non-redundant annotation
    res_mat_mod <- res_mat_filt[rownames(res_mat_filt) %in% tr2g$gene,]
    
    # Change gene IDs to gene symnbols
    rownames(res_mat_mod) <- tr2g$gene_name[match(rownames(res_mat_mod), tr2g$gene)]
    
    # Create Seurat object
    print("Creating Seurat object...")
    seu_tmp <- CreateSeuratObject(res_mat_mod, min.cells = min.cells, min.features = min.features, project = samp_nm)
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
      print("Filter out doublets...")
      seu_tmp <- subset(seu_tmp, subset = DubFind == "Singlet")
    }
    
    
    
    # Get cell numbers for data.frame
    og_cells <- dim(res_mat_filt)[2]
    filt_cells <- dim(seu_tmp)[2]
    
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


