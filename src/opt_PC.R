##########################
# This function determins the optimal PC for PCA using the elbow method


##########################
# https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
# Optimal PC is the lesser of the two test below
# The point where the principal components only contribute 5% of standard deviation and the principal components cumulatively contribute 90% of the standard deviation.
# The point where the percent change in variation between the consecutive PCs is less than 0.1%.
#########################



opt_PC <- function(seu_obj, assay = "RNA"){
  
  
  ##########################
  # This was written in Sept 29th, 2022 by MJB
  # USAGE: opt_PC(seu_obj, assay = "RNA")
  #
  # seu <- Seurat object
  # assay <- assay for using in PCA, default = "RNA"
  ##########################
  
  
  # Plot the elbow plot
  ElbowPlot(object = seu_obj, 
            ndims = 50)
  
  # Determine percent of variation associated with each PC
  pct <- seu_obj[["pca"]]@stdev / sum(seu_obj[["pca"]]@stdev) * 100
  
  # Calculate cumulative percents for each PC
  cumu <- cumsum(pct)
  
  # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
  co1 <- which(cumu > 90 & pct < 5)[1]
  
  co1
  
  # Determine the difference between variation of PC and subsequent PC
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  
  # last point where change of % of variation is more than 0.1%.
  co2
  
  # Minimum of the two calculation
  pcs <- min(co1, co2)
  
  pcs
  
  # Create a dataframe with values
  plot_df <- data.frame(pct = pct, 
                        cumu = cumu, 
                        rank = 1:length(pct))
  
  # Elbow plot to visualize 
  ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
    geom_text() + 
    geom_vline(xintercept = 90, color = "grey") + 
    geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
    theme_bw()
  
  return(pcs)
  
}