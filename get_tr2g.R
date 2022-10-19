##########################
# This is a helper function for getting transcript to gene information from a gtf file


get_tr2g <- function(gtf_file){
  
  
  ##########################
  # This was written in Sept 29th, 2022 by MJB
  # USAGE: tr2g(gft_file = "path_to_gtf")
  #
  # gtf_file <- path to the gtf file
  ##########################
  
  
  
  require(GenomicFeatures)
  
  # Import gtf as data.frame and get the base Gene ID name
  gtf.gr <- rtracklayer::import(gtf_file)
  gtf.df <- as.data.frame(gtf.gr)
  gtf.df$gene_id <- gsub("\\..+", "", gtf.df$gene_id)
  
  # Generate the transcript to gene data.frame
  tr2g <- unique(gtf.df[,c("gene_id", "gene_name")])
  tr2g <- tr2g[order(tr2g$gene_id),]
  tr2g <- tr2g[!duplicated(tr2g$gene_name),]
  colnames(tr2g) <- c("gene", "gene_name")
  
  return(tr2g)
  
}




# TESTING ####################
# tr2g(gtf_file = "../data/external/gencode.vM25.annotation.gtf")


