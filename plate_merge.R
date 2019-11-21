library(dplyr)

#' Title
#'
#' @return
#' @export
#'
#' @examples
plate_merge <- function(outdir, id_list){
  id_list <- id_list %>% strsplit(', ') %>% unlist %>% 
    gsub("^.*_(.*)\\..*$","\\1",.) %>% paste0("_matrix.tsv")
  
  quant_dir_content <- list.files(outdir, full.names = T)
  file_list <- sapply(id_list, function(p) grep(p, quant_dir_content, value = T)) 
  class(file_list)
  
  tables_list <- lapply(file_list,
         FUN=function(files){read.table(files, check.names = F,
                                        header=T, sep="\t")
         }
  )
  names(tables_list) <- NULL
  
  final_matrix <- do.call("cbind", tables_list)
  
  write.table(final_matrix, paste0(outdir, "all_plates.tsv"), sep = '\t')
}


if(sys.nframe()==0){
  outloc <- "./"
  plate_info <- 0
  sheetfile <- NULL
  mtrds <- NULL
  
  args=commandArgs(trailingOnly=T)
  
  if(!file.exists(args[1]))
    stop(paste0("Matrix file ('",args[1],"') not found!"))
  
  if(length(args)>=2){
    outloc <- args[1]
    ids <- args[2]
    if(!file.exists(outloc)){
      warning(paste0("Output directory ('",outloc,"') does not exist!"),
              "\nCreating the directory...")
      system(paste("mkdir -p",outloc))
    }
  } else {
    stop(">>> Missing args! <<<")
  }
  
  plate_merge(args[1], args[2])
}