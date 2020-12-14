# this script requires internet connection

library(tidyverse)
library(biomaRt)

# human:"hsapiens_gene_ensembl"
# mouse:"mmusculus_gene_ensembl"

function(species="hsapiens_gene_ensembl",
         biomart_attributes=c("ensembl_gene_id","ensembl_transcript_id_version",
                              "external_transcript_name", "external_gene_name")){
  
  mart_obj <- useMart("ensembl", dataset = mart_species)
  db <- getBM(biomart_attributes, mart=mart_obj)
  saveRDS(paste(species, "mart","db.rds", sep = "_"), db)
}

