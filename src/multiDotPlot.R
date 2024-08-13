multiDotPlot <- function(dat_seu, feat, xval, xval_ord, grp, grp_ord, assay_exp="SCT", cnt_exp="data", assay_pop="RNA", cnt_pop="data", col_pal=NULL, title_name){
  
  #############################################
  # Written by MJB on Aug 13th, 2024, last modified on 
  # dat_seu = seurat object
  # feat = vector of features for dotplot (genes)
  # xval = seu@meta.data column name for the dotplot x-axis
  # xval_ord = order of features on x-axis, character vector
  # grp = seu@meta.data column name for the group name (facet)
  # grp_ord = order for the group name (facet), character vector
  # assay_exp = assay for expression data
  # cnt_exp = slot in assay from which to pull expression data
  # assay_pop = assay for population data
  # cnt_pop = slot in assay from which to pull population data
  # col_pal = color palette to use for fill
  # title_name
  #
  # Code modified from https://davemcg.github.io/post/lets-plot-scrna-dotplots/ 
  #############################################

  library(Seurat)
  library(tidyverse)

  if (is.null(col_pal)){
    col_pal <- rev(RColorBrewer::brewer.pal(11, "RdYlBu"))
  }
  
  # Initialize expression matrix
  tmp_exp_mat <- data.frame(row.names = feat)        
  
  # Get data.frame of average expression for each age/subcell type 
  for (ii in levels(dat_seu@meta.data[[xval]])){
    
    # Get average expression
    
    expr <- FetchData(object = dat_seu, vars = xval)
    tmp_seu <- dat_seu[, which(x = expr == ii)]
    
    # tmp_seu <- subset(dat_seu, subset = `xval` == ii)
    tmp_mat <- AverageExpression(tmp_seu, group.by = grp,
                                 features = feat,
                                 assays = assay_exp, slot = cnt_exp)
    
    # Tidy the DF
    tmp_mat_tidy <- tmp_mat[[1]] %>% 
      as.data.frame() %>%
      rownames_to_column(var = "Gene") %>% 
      gather("CellType", "Count", -Gene)
    
    # Add ages to cell type
    tmp_mat_tidy <- tmp_mat_tidy %>% 
      mutate(CellType_age = paste(CellType, ii, sep = "_"),
             Age = ii)
    
    
    # Get number of cells expressing the genes per SubCellType at this age
    tmp_cell_exp_ct_df <- data.frame()
    for (iii in unique(tmp_mat_tidy$CellType)){
      
      expr <- FetchData(object = tmp_seu, vars = grp)
      tmp_celltype_boul <- (expr == iii) %>% as.vector()
      
      tmp_dat <- GetAssayData(tmp_seu, assay = assay_pop, slot = cnt_pop)[feat,tmp_celltype_boul]
      tmp_sum <- tmp_dat %>% as.data.frame() %>% apply(1, function(x){sum(x>0)})
      tmp_cell_exp_ct_df <- rbind(tmp_cell_exp_ct_df,
                                  data.frame(cellType = iii,
                                             cell_ct = dim(tmp_dat)[2],
                                             cell_exp_ct = tmp_sum))
    }
    
    tmp_mat_tidy <- tmp_mat_tidy %>% 
      mutate(cell_exp_ct = tmp_cell_exp_ct_df$cell_exp_ct,
             cell_ct = tmp_cell_exp_ct_df$cell_ct)
    
    # rbind the different timepoints
    tmp_exp_mat <- rbind(tmp_exp_mat, tmp_mat_tidy)
    # ord <- paste(tmp_exp_mat$CellType %>% unique() %>% rep(each = 4), seu_muel$age %>% unique(), sep = "_")
    # tmp_exp_mat$CellType_age <- factor(tmp_exp_mat$CellType_age, levels = ord)
  } 

  
  # Add gene-wise Z-score
  tmp_exp_mat <- tmp_exp_mat %>% 
    group_by(Gene) %>% 
    mutate(Zscore = scale(Count, scale = T) %>% unlist() %>% as.vector())
  
  # Add plot data to list
  df <- tmp_exp_mat
  
  df$Gene <- factor(df$Gene, levels = rev(df$Gene %>% unique()))
  
  
  # Generate plot
  
  tmp_plot <- df %>% 
    
    mutate(`% Expressing` = (cell_exp_ct/cell_ct) * 100,
           Age = factor(Age, levels = xval_ord),
           Gene = factor(Gene, levels = rev(feat)),
           Grp = factor(CellType, levels = grp_ord)) %>% 
    filter(Count > 0, `% Expressing` > 1) %>% 
    
    ggplot(aes(x=Age, y = Gene, group = Grp, color = Zscore, size = `% Expressing`)) + 
    geom_point() +
    scale_size_continuous(limits = c(0,100)) +
    theme(axis.line  = element_blank()) +
    theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1)) +
    ylab('') + xlab('') +
    theme(panel.grid.major = element_blank(),
          panel.background = element_rect(fill = "grey95")) +
    scale_color_gradientn(colours = col_pal,
                          # scale_color_gradientn(colours = rev(vir_pal),
                          # limits = c(-3,3), oob = scales::squish,
                          name = 'Average Expression') + 
    ggtitle(title_name) +
    facet_grid(.~Grp, scales = "fixed") 
  
  return(tmp_plot)
}


# Testing
# seu.rods <- readRDS(file = "/Volumes/Accelsior_8M2/scRNAseq/scRNA_aging/cellranger/data/interim/seu_rod_final_230808.rds")
# seu.rods$seurat_clusters <- seu.rods$SCT_snn_res.0.8
# seu.rods$seurat_clusters <- factor(seu.rods$seurat_clusters, levels = as.character(c(20,14,8,4,33,3,5,1,0)))
# rod_clst <- ifelse(seu.rods$seurat_clusters %in% c("4", "8", "14", "20"), "visTrans", "synEnr")
# seu.rods <- AddMetaData(object = seu.rods, metadata = rod_clst, col.name = "rod_cluster")
# seu.rods$rod_cluster <- factor(seu.rods$rod_cluster, levels = c("visTrans", "synEnr"))
# dat_seu <- seu.rods
# 
# feat <- c("Rho", "Nrl", "Gnat1", "Pde6b", "Sag", "Rcvrn", "Rpgrip1", "Syt1")
# xval <- "age"
# xval_ord <- dat_seu$age %>% unique()
# grp <- "rod_cluster"
# grp_ord <- dat_seu$rod_cluster %>% unique() %>% rev()
# assay_exp="SCT"
# cnt_exp="data"
# assay_pop="RNA"
# cnt_pop="data"
# col_pal=NULL
# title_name=NULL


