# if run as part of the QC pipeline
# output directory is set within NF process (to quantifications dir)

library(biomaRt)
library(tidyverse)


tx2gene_counts <- function(counts_table_location, outloc, species="Hsapiens",
                           hdb="../References/human_mart_db.rds",
                           mdb="../References/mouse_mart_db.rds"){
  
  c_matrix <- read.table(counts_table_location, header = T, sep='\t',
                         check.names = F)
  
  if(species=="Hsapiens"){
    mart_species <- "hsapiens_gene_ensembl"
    db <- readRDS(hdb)
  } else if(species=="Mmusculus"){
    mart_species <- "mmusculus_gene_ensembl"
    db <- readRDS(mdb)
  } else {
    print(species)
    stop("<--- Invalid species string! --->")
  }
  
  # mart_obj <- useMart("ensembl", dataset = mart_species)
  # db <- getBM(c("ensembl_gene_id","ensembl_transcript_id_version",
  #               "external_transcript_name", "external_gene_name"), mart=mart_obj)
  
  # matching rownames to biomart gene names
  ids <- rownames(c_matrix) %>% gsub("\\|.*$","",.)
  c_matrix$gene_name <- db[match(ids, db$ensembl_transcript_id_version),]$external_gene_name
  
  # removing transcripts without an associated gene name
  c_matrix <- c_matrix[-which(is.na(c_matrix$gene_name)),]
  
  # R magic
  gene_level_matrix <- aggregate(. ~ gene_name, c_matrix, sum)
  rownames(gene_level_matrix) <- gene_level_matrix$gene_name
  
  write.table(gene_level_matrix[,-1], paste0(outloc,"plates_as_genelevel.tsv"), sep='\t')
}



if(sys.nframe()==0){

  args=commandArgs(trailingOnly=T)
  
  
  if(!file.exists(args[1]))
    stop(paste0("<--- Matrix file ('",args[1],"') not found! --->"))
  
  if(!file.exists(args[2]))
    stop(paste0("<--- Tx2g output location ('",args[2],"') not found! --->"))
  
  if(length(args)<2)
    stop("<--- Missing arguments in tx2g call! --->")
  
  if(length(args)>=3){
    species <- args[3]
    tx2gene_counts(args[1], args[2], species)
  } else {
    tx2gene_counts(args[1], args[2])
  }
  
  
}

