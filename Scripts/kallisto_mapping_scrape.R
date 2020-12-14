library(dplyr)
library(rjson)

#' Extracts pseudoalignment percentages as a table of sample names and a mean
#' from a folder with all the sample quantifications subfolders
#'
#' @param quant_superfolder Name of the main folder. Typically intermediate output of the QC pipeline.
#' @param output_name Name of the file containing alignment percentages and their mean.
#' Defaults to 'percent_pseudoaligned.txt'
#' @param metric Defaults to 'p_pseudoaligned'. If not found in 'run_info.json', it will throw an error/
#' @param qc_name Defaults to '/qc'. Name that needs to be ignored in the folder.
#' If result of the QC pipeline, it's probably there.
#'
#' @return Nothing. Writes a small file.
#' @export
#'
#' @examples 
kallisto_mapping_scrape <- function(quant_superfolder, output_name="percent_pseudoaligned.txt",
                                    metric="p_pseudoaligned", qc_name="/qc"){
  paths <- list.files(quant_superfolder, full.names = T)
  paths <- paths %>% dir.exists %>% paths[.]
  qc_fold <- grep(qc_name, paths)
  if(length(qc_fold)>0)
    paths <- paths[-qc_fold]
  json_files <- paste0(paths, "/run_info.json")
  
  sapply(json_files, function(jf) paste(readLines(jf), collapse = "") %>%
           fromJSON() %>% .[[metric]]) -> p_align_per_well
  
  names(p_align_per_well) %>%
    substring(1,nchar(.)-nchar("/run_info.json")) %>%
    gsub("^.*/","",.) -> names(p_align_per_well)
  
  p_align_per_well <- c("mean"=mean(p_align_per_well),p_align_per_well)
  
  write.table(p_align_per_well %>% as.data.frame, output_name, col.names = F)
}


if(sys.nframe()==0){
  args=commandArgs(trailingOnly=T)
  if(!file.exists(args[1]))
    stop(">>> Folder location not found! <<<")
  if(length(args)<2)
    stop(">>> Output name not provided! <<<")
  runinfo_metric <- ifelse(length(args)>2, args[3],"p_pseudoaligned")
  
  kallisto_mapping_scrape(quant_superfolder=args[1], output_name=args[2], metric=runinfo_metric)
}