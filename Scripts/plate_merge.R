library(dplyr)

#' Title
#'
#' @return
#' @export
#'
#' @examples
plate_merge <- function(outdir, q_merge_col='est_counts', id_list){
  # v1: id_list
  # > id_list
  # [1] "[/ei/cb/development/lany/CB-Iain_Macaulay_EI_CU_ENQ-5072_A_01/Analysis/scqc_reqs-0.2.3/test2.4plate/work/fa/2b420fbaa97dd8af0244eb9acca1ba/Finished_CU5DAY0.txt, /ei/cb/development/lany/CB-Iain_Macaulay_EI_CU_ENQ-5072_A_01/Analysis/scqc_reqs-0.2.3/test2.4plate/work/9d/2fe712a80b3945ca6c51cd0ded05e4/Finished_CU7DAY0.txt, /ei/cb/development/lany/CB-Iain_Macaulay_EI_CU_ENQ-5072_A_01/Analysis/scqc_reqs-0.2.3/test2.4plate/work/87/4592ac35208e669b988678ded0f1b7/Finished_CU5DAY7.txt, /ei/cb/development/lany/CB-Iain_Macaulay_EI_CU_ENQ-5072_A_01/Analysis/scqc_reqs-0.2.3/test2.4plate/work/b8/3313f5215f177dcf91ecd90816882f/Finished_CU7DAY7.txt]"
  # 
  id_list <- id_list %>% strsplit(' ') %>% unlist %>% gsub('\\[','',.) %>% gsub('\\]','',.) %>% basename() %>%
    sub('\\..*$', '',.) %>% gsub('^.*_','',.) %>% paste0(q_merge_col,.,"_matrix.tsv")
  file_list <- file.path(outdir,id_list)
  
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
  
  if(length(args)==3){
    outloc <- args[1]
    if(!file.exists(outloc)){
      warning(paste0("Output directory ('",outloc,"') does not exist!"),
              "\nCreating the directory...")
      system(paste("mkdir -p",outloc))
    }
    plate_merge(args[1], args[2], args[3])
  } else {
    stop(">>> Missing args! <<<")
  }
}