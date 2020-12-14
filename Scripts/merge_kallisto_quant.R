
#' Merges abundance files from kallisto quantifications
#'
#' @param folders_loc Folder containing folders, each
#' being a result of one `kallisto quant`
#' @param metric Which `abundance.tsv` column to extract.
#' Default is `tpm`. Other options - `est_counts`, `length`,
#' `eff_length` and `target_id`.
#'
#' @return Nothing. Writes a merged tsv in  folders_loc.
#' @export
#'
#' @examples
merge_kallisto_quant <- function(folders_loc="~/AThaliana/kallisto_q/",
                                 metric="tpm", plate_id=NULL){

	kquant_folds <- list.dirs(folders_loc, full.names = F)
	kquant_folds <- unique(gsub("/.*$","",kquant_folds))
	qcfold_ind <- grep("qc",kquant_folds)	
	if(length(qcfold_ind)>0){
		kquant_folds <- kquant_folds[-qcfold_ind]
	}
	kquant_folds <- kquant_folds[kquant_folds!=""]
	if(!is.null(plate_id)){
	  # filtering folders for the specified plate
	  kquant_folds <- grep(plate_id, kquant_folds, value = T)
	}
	sample_names <- kquant_folds
	abundance_paths <- paste0(folders_loc,kquant_folds, "/abundance.tsv")

	singleton <- read.table(abundance_paths[1], header=T)
	index <- which(colnames(singleton)==metric)

	if(identical(index, integer(0)))
	    stop(paste("No column named",metric))

	tpm_matrix <- do.call("cbind",
		           lapply(abundance_paths,
		                          FUN=function(files){read.table(files,
		                                                         header=TRUE, sep="\t")[index]
		                            }
		                  )
	)

	colnames(tpm_matrix) <- sample_names
	rownames(tpm_matrix) <- singleton$target_id
	
	
	# adding plate label to the counts file for selection/distinction later
	if(is.null(plate_id)){
		plate_label <- ""
	}
	else{
		plate_label <-  as.character(plate_id)
	}

	write.table(tpm_matrix, paste0(folders_loc,metric,
	                               plate_label,"_matrix.tsv"), sep = "\t")
}

if(sys.nframe()==0){
	args=commandArgs(trailingOnly=T)
	if(!file.exists(args[1]))
		stop(">>> Folders location not found! <<<")
	column_metric <- ifelse(length(args)>1, args[2],"tpm")
  pid <- ifelse(length(args)>2, args[3], NULL)
	
	merge_kallisto_quant(args[1], column_metric, pid)
}



#' Title
#'
#' @param tag 
#'
#' @return
#' @export
#'
#' @examples
revcom <- function(tag){
  
  compl <- function(dna_b){
    case_when(
      dna_b == 'A' ~ 'T',
      dna_b == 'C' ~ 'G',
      dna_b == 'G' ~ 'C',
      dna_b == 'T' ~ 'A',
      TRUE ~ dna_b
    )
  }
  
  if(nchar(tag)>1){
  tag %>% strsplit('') %>% unlist() %>% compl %>% rev %>% paste0(collapse = '')
  }else{
    tag
  }
}
