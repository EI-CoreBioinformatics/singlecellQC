# this script requires internet connection

library(tidyverse)
library(biomaRt)

# human:"hsapiens_gene_ensembl"
# mouse:"mmusculus_gene_ensembl"

make_biomart_db <- function(mart_species="hsapiens_gene_ensembl",
         biomart_attributes=c("ensembl_gene_id","ensembl_transcript_id_version",
                              "external_transcript_name", "external_gene_name")){
  
  mart_obj <- useMart("ensembl", dataset = mart_species)
  db <- getBM(biomart_attributes, mart=mart_obj)
  saveRDS(db,paste(mart_species, "mart","db.rds", sep = "_"))
}

if(sys.nframe()==0){
  args=commandArgs(trailingOnly=T)
  
  if(length(args)<1)
    stop("No specified species!")
  if(length(args)<2){
    make_biomart_db(args[1])
   } else{
     make_biomart_db(args[1], args[2])
  }
}