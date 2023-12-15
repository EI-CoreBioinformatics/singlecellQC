#library(biomaRt)
# library(tidyverse)
library(dplyr)

tx2gene_counts <- function(counts_table_location, outloc, species="Bos_taurus",
                           dbfile="/ei/cb/development/lany/CB-Iain_Macaulay_EI_CU_ENQ-5072_A_01/Scripts/scqc/References/btaurus_gene_ensembl_mart_db.rds")
{
  c_matrix <- read.table(counts_table_location, header = T, sep='\t',
                         check.names = F)
  
  if(species=="Hsapiens" || species=="Bos_taurus"){
    # mart_species <- "hsapiens_gene_ensembl"
    db <- readRDS(dbfile)
  } else {
    mapping_matrix <- read.table(dbfile, header=F, sep='\t', check.names=F,stringsAsFactors = F)
    colnames(mapping_matrix) <- c('transcript', 'gene')
  }
  
  # matching rownames to biomart gene names
  if (species=="Hsapiens") {
    # rowname: ENST00000456328.2|ENSG00000223972.5|OTTHUMG00000000961.2|OTTHUMT00000362751.1|DDX11L1-202|DDX11L1|1657|processed_transcript|
    ids <- rownames(c_matrix) %>% gsub("\\|.*$","",.)
    c_matrix$gene_name <- db[match(ids, db$ensembl_transcript_id_version),]$external_gene_name
    # removing transcripts without an associated gene name
    c_matrix <- c_matrix[-which(is.na(c_matrix$gene_name)),]
  } else if (species=="Bos_taurus"){
    # rowname: "transcript:ENSBTAT00000007786" ; ids: ENSBTAT00000007786
    ids <- rownames(c_matrix) %>% gsub("transcript:","",.)
    # head(db.cow$ensembl_transcript_id_version)
    # [1] "ENSBTAT00000008737.6" wanna remove ".6"
    tmp.ensembl_transcript_id_version <- str_replace_all(db$ensembl_transcript_id_version,
                                                        "\\.\\d+",replacement="")
    # c_matrix$gene_name <- db.cow[match(ids, db.cow$ensembl_transcript_id_version),]$external_gene_name
    c_matrix$gene_name <- db[match(ids, tmp.ensembl_transcript_id_version),]$external_gene_name
    # removing transcripts without an associated gene name
    c_matrix <- c_matrix[-which(c_matrix$gene_name==""),]
  } else {
    # c_matrix rowname: "transcript:ENSBTAT00000007786" ; ids: ENSBTAT00000007786
    ids <- rownames(c_matrix) %>% gsub("transcript:","",.)
    c_matrix$gene_name <- mapping_matrix[match(ids, mapping_matrix$transcript),]$gene 
    # removing transcripts without an associated gene name
    if (length(which(c_matrix$gene_name==""))>0){
      c_matrix <- c_matrix[-which(c_matrix$gene_name==""),]
    }
  }
  
  # R magic
  gene_level_matrix <- aggregate(. ~ gene_name, c_matrix, sum)
  rownames(gene_level_matrix) <- gene_level_matrix$gene_name
  
  write.table(gene_level_matrix[,-1], paste0(outloc,"plates_as_genelevel.tsv"), sep='\t')
}



if(sys.nframe()==0){
  # est_counts_tx2gene.lany.R ./quants_dir/est_countsCU5DAY0_matrix.tsv ./quants_dir/ Bos_taurus \
  #    ./References/btaurus_gene_ensembl_mart_db.rds
  args=commandArgs(trailingOnly=T)
  
  
  if(!file.exists(args[1]))
    stop(paste0("<--- Matrix file ('",args[1],"') not found! --->"))
  
  if(!file.exists(args[2]))
    stop(paste0("<--- Tx2g output location ('",args[2],"') not found! --->"))
  
  if(length(args)<2)
    stop("<--- Missing arguments in tx2g call! --->")

  if(length(args)<3)
    stop("<--- Missing species tx2g call! --->")

  # if(length(args)>=3){
  #   species <- args[3]
  #   tx2gene_counts(args[1], args[2], species)
  # } else {
  #   tx2gene_counts(args[1], args[2])
  # }

  if(length(args)>=4){
    species <- args[3]
    dbfile <- args[4]
    tx2gene_counts(args[1], args[2], species, dbfile)
  } else {
    tx2gene_counts(args[1], args[2], species)
  }
  
}
