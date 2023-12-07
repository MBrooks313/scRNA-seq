##########################
# This is a helper function for importing kallisto scRNA-seq files into a matrix


read_count_output <- function(dir, name="cells_x_genes") {
  
  
  ##########################
  # This was written in Sept 29th, 2022 by MJB
  # USAGE: read_count_output(dir, name = "cells_x_genes")
  #
  # dir <- path to the sample kallisto output directory
  # name <- result matrix name to import 
  ##########################
  
  
  
  require(Matrix)
  
  dir <- normalizePath(dir, mustWork = TRUE)
  m <- readMM(paste0(dir, "/", name, ".mtx"))
  m <- Matrix::t(m)
  m <- as(m, "dgCMatrix")
  # The matrix read has cells in rows
  ge <- ".genes.txt"
  genes <- readLines(file(paste0(dir, "/", name, ge)))
  barcodes <- readLines(file(paste0(dir, "/", name, ".barcodes.txt")))
  colnames(m) <- barcodes
  rownames(m) <- genes
  return(m)
}