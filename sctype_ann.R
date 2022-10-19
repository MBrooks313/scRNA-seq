##########################
# This function performs annotation on Seurat object


sctype_ann <- function(seu, gs_pos=NULL, gs_neg=NULL, sctissue="Eye", sctypes=NULL, assay_data="RNA"){
  
  
  ##########################
  # This was written in Sept 29th, 2022 by MJB
  # USAGE: sctype_ann(seu, gp_pos=<cell-goi_list>, gs_neg, sctissue="Eye", sctypes=NULL, assay_data="RNA")
  #
  # seu <- Seurat object for annotating
  # gs_pos <- list of cell types with genes for positive association 
  # gs_neg <- list of cell types with genes for negative association
  # tissue <- tissue to pull scType db cell annotations 
  # sctypes <- additional cell types pulled from scType db 
  # assay_scale <- assay scaled.data slot used for cell type scoring, c("RNA", "SCT", "integrated") 
  #
  # Value
  # seu <- Seurat object with annotation
  ##########################


  require("dplyr")
  require("Seurat")
  require("HGNChelper")
  
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
  
  
  # DB file
  db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
  tissue = sctissue # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain
  
  # prepare scType db gene sets
  gs_list = gene_sets_prepare(db_, tissue)
  
  # Use sc_type db if none provided
  if (is.null(gs_pos)){
    gs_pos <- gs_list$gs_positive
  }
  
  # Add scTpye db cell types if needed
  if (!is.null(sctypes)){
    for (i in 1:length(sctypes)){
      gs_pos[[sctypes[i]]] <- gs_list$gs_positive[[sctypes[i]]]
    }
  }
  
  
  # get cell-type by cell matrix
  es.max = sctype_score(scRNAseqData = seu[[assay_data]]@scale.data, scaled = TRUE, 
                        gs = gs_pos, gs2 = gs_neg) 
  
  # merge by cluster
  cL_resutls = do.call("rbind", lapply(unique(seu@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(seu@meta.data[seu@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seu@meta.data$seurat_clusters==cl)), 10)
  }))
  sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
  
  # set low-confident (low ScType score) clusters to "unknown"
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
  print(sctype_scores[,1:3])
  
  seu@meta.data$CellType = ""
  for(j in unique(sctype_scores$cluster)){
    cl_type = sctype_scores[sctype_scores$cluster==j,]; 
    seu@meta.data$CellType[seu@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
  }
  
  return(seu)

}
