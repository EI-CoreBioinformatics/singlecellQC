
library(dplyr)
library(scater)

#' Running single-cell QC on a tpm/count matrix
#'
#' @param tsv_location Expression matrix
#' @param output_location Where to store the QC files.
#' @param plate_info 0 for a single plate like the first GP testrun,
#' 1 for Anita/Geoff multi-plate,
#' 2 (default) for sample_sheet provided ones
#' @param sheet_file Sample sheet with plate position, control status and sample ID information
#' If missing, plate position will try to be inferred from the .fastq file names
#' @param mt_file Mitochondrial genes. If not provided, the ones starting with "MT-" or
#' "mt-" will be considered mitochondrial (works only for human and mouse)
#'
#' @return None. Generates files (pdf plots, tsv tables).
#' @export
#'
#' @examples
scqc_from_tsv <- function(tsv_location, output_location="./", plate_info=2,
                          sheet_file=NULL, mt_file=NULL) {
  #tsv will be est_counts or tinds_cellctrlpm matrix transcript X sample
  #but counts still used as label in th SCE object
  counts_df <- read.table(tsv_location, header = T, check.names = F)
  
  experiment_id <- gsub("/[A-z|_|\\.]+$","",tsv_location) %>% 
    gsub("^.*/","",.)
  plate_id <- "1"
  
  # using the sheet_file if provided
  if(!is.null(sheet_file)){
    samplesheet <- read.csv(sheet_file)
    # samplesheet column names assumed standardised

    row_ids <- sapply(samplesheet$unique_sample_id_suffix,
           function(id) grepl(id, colnames(counts_df))) %>% apply(1,which)
    
    well_id <- samplesheet[row_ids,]$well
    cell_number <- samplesheet[row_ids,]$number_of_cells
    inds_cellctrl <- which(samplesheet[row_ids,]$control==T)
    
    # could this not be unique?
    # one plate, multiple experiments?
    experiment_id <- samplesheet[row_ids,]$experiment %>% unique
    
    plate_id <- samplesheet[row_ids,]$plate_id %>% unique
    
    colnames(counts_df) <- samplesheet[row_ids,]$unique_sample_id_suffix
    
    # cell_type and cell_phenotype probably best not too always excpect?
  }
  #print(plate_id)
  if(length(plate_id)>1){
    plate_id <- "all"
  }
  
  #feature and cell controls
  inds_ercc <- counts_df %>% rownames() %>% grep("^ERCC",.)
  inds_mito <- counts_df %>% rownames() %>% grep("\\|MT-",., ignore.case = T) #works only on HS
  if(!exists("inds_cellctrl")) inds_cellctrl <- numeric(0)
  
  long_names <- counts_df[-inds_ercc,] %>% rownames()
  #takes a while
  short_names <- sapply(long_names, function(s) strsplit(s,"\\|") %>% .[[1]] %>% .[1]) %>% as.vector 
  rownames(counts_df)[-inds_ercc] <- short_names

  mito_names <- long_names[inds_mito] %>% strsplit("\\|") %>% sapply(function(x) x[5]) %>% as.vector
  rownames(counts_df)[inds_mito] <- mito_names
  
  if(!is.null(mt_file)){
    mito_names <- readRDS(mt_file)
    inds_mito <- counts_df %>% rownames %>%
      match(mito_names, .) %>% na.omit()
  }
    
  
  sce <- SingleCellExperiment(assays=list(counts=counts_df %>% as.matrix))
  sce <- calculateQCMetrics(sce,
                            feature_controls = list(ERCC = inds_ercc, mitochondrial = inds_mito),
                            cell_controls = list(Controls = inds_cellctrl) 
                              )
  
  rm(counts_df)
  setwd(output_location)
  

  
  
  #total count across genes
  plot_highest_expression <- plotHighestExprs(sce, exprs_values = "counts")
  #mean expression vs frequncy expressed
  plot_expr_vs_mean <-  plotExprsFreqVsMean(sce, exprs_values = "counts", show_smooth = F)
  #cumulative distribution plot
  
  cellctrl_flag <- NULL
  if(length(inds_cellctrl)>0){
    cellctrl_flag <- "is_cell_control"
    }
  
  plot_cumulative_dist <- plotScater(sce, nfeatures = 300, exprs_values = "counts",
                                     colour_by = cellctrl_flag)  

  #%MT - 
  plot_mito_scatter <- plotColData(sce, x = "total_features_by_counts",
              y = "pct_counts_mitochondrial",
              colour_by = cellctrl_flag) +
    theme(legend.position = "top") + xlab("total_features") +
    labs(title="Mitochondrial expression percentage ~ total features")
  #%ERCC
  plot_ercc_scatter <- plotColData(sce, x = "total_features_by_counts",
              y = "pct_counts_ERCC",
              colour_by = cellctrl_flag) +
    theme(legend.position = "top") + xlab("total_features") +
    labs(title="Spike-in expression percentage ~ total features")
  #counts across cells
  #% any feature control
  plot_fctrl_scatter <- plotColData(sce, x = "total_features_by_counts",
              y = "pct_counts_feature_control",
              colour_by = cellctrl_flag) +
    xlab("total_features") +
    #stat_smooth(method = "lm", se = FALSE, size = 1.5, fullrange = TRUE) +
    ggtitle(paste("Out of", sce@assays$data$counts %>% nrow, "transcripts"))
  
  #counts ~ total
  plot_featurenum_vs_counts <- plotColData(sce, x = "log10_total_counts",
              y = "total_features_by_counts",
               colour_by = cellctrl_flag) +
    theme(legend.position = "top") + xlab("log_total_counts") + ylab("total_features") +
    labs(title="Total features ~ total counts (log scale)")
  
  pdf(paste0("QC_scatterplots_",plate_id,".pdf"), pointsize = 12)
  plot_featurenum_vs_counts %>% plot
  plot_cumulative_dist %>% plot
  plot_expr_vs_mean %>% plot
  plot_highest_expression %>% plot
  
  plot(plot_fctrl_scatter)
  plot(plot_mito_scatter)
  plot(plot_ercc_scatter)
  dev.off()
  
  
  #0 - 1st 2 GP test runs naming
  #1 - Anita/Geoff naming
  #2 - run 3 (Jim's samplesheet) GP naming
  sce$Cell <- sce@assays$data$counts %>% colnames
  if(exists("well_id")) sce$plate_position <- well_id
  if(plate_info==0){
  sce$plate_position <- sce$Cell %>% gsub("^([[:alnum:]]+)\\_.*$","\\1",.) %>%
    {ifelse(nchar(.)==2,gsub("(.)(.)","\\10\\2",.),.)}
  } else if (plate_info==1){
    c_string <- sce$Cell
    sce$plate_position <- c_string %>% regexec("plate",., ignore.case = T) %>%
      .[[1]] %>% .[1] %>% `+`(nchar("plateX")) %>% 
      substring(c_string, .) %>%  #same format as standard 1 plate GP now
      gsub("^([[:alnum:]]+)\\_.*$","\\1",.) %>%
      {ifelse(nchar(.)==2,gsub("(.)(.)","\\10\\2",.),.)}
      
  } else if(plate_info==2){
    sce$plate_position <- well_id
  } else {
    stop("plate_info error")
  }
  
  plot_plate_mito <- plotPlatePosition(sce,
                    colour_by = "pct_counts_mitochondrial",
                    by_exprs_values = "counts",
                    point_size = 9) +
    guides(fill="colorbar") + theme(legend.position = "top")
  plot_plate_ercc <- plotPlatePosition(sce,
                                       colour_by = "pct_counts_ERCC",
                                       by_exprs_values = "counts",
                                       point_size = 9) +
    guides(fill="colorbar") + theme(legend.position = "top")
  plot_plate_counts <- plotPlatePosition(sce,
                                         colour_by = "log10_total_counts",
                                         by_exprs_values = "counts",
                                         point_size = 9) + guides(fill="colorbar")
  # counts per well
  pdf(paste0("Plate_position_plots_",plate_id,".pdf"), pointsize = 12)
  plot(plot_plate_counts)
  plot(plot_plate_ercc)
  plot(plot_plate_mito)
  dev.off()
  
  # ERCC mol
  # >add %ERCC
  ercc_counts <- sce@assays$data$counts[inds_ercc,] %>% as.data.frame(stringsAsFactors=F)
  colnames(ercc_counts) <- ercc_counts %>% colnames() %>% gsub("_.*$","",.)
  ercc_counts$All_samples <- rowSums(ercc_counts)
  ercc_counts$ERCC <- rownames(ercc_counts)
  df_col <- ncol(ercc_counts)
  
  ercc_table <- ercc_counts[,c(df_col, df_col-1, 1:(df_col-2))]
  
  write.table(ercc_table, paste0("Sample_counts_across_ERCCs_",plate_id,".tsv"),
              sep='\t', row.names = F)

   # sce$cell_type <- sce %>% colData() %>% rownames %>%
   #   match(samplesheet$unique_sample_id_suffix) %>%
   #   samplesheet[.,] %>% `$`(cell_type)

  sce <- normalize(sce)
  sce <- runPCA(sce, ncomponents=20, detect_outliers=T)
  sce <- runTSNE(sce)
  
  control_shape <- NULL
  if(length(inds_cellctrl)) control_shape <- "is_cell_control"
  
  control_scale <- NULL
  if(length(inds_cellctrl)) control_scale <- scale_shape_manual(values=c(16,9))
  
  plot_pca <- plotPCA(sce,
                      colour_by = "log10_total_counts",
                      size_by = "total_features_by_counts",
                      shape_by=control_shape) + control_scale
  
  plot_tsne <- plotTSNE(sce,
           colour_by = "pct_counts_feature_control",
           size_by = "total_features_by_counts",
           shape_by = control_shape) + control_scale
  
  # plot_tsne_type <- plotTSNE(sce,
  #                       colour_by = "cell_type",
  #                       size_by = "total_features_by_counts",
  #                       shape_by = "cell_type")
  
  pdf(paste0("QC_dimReductions_",plate_id,".pdf"), pointsize=12)
  plot(plot_pca)
  plot(plot_tsne)
  dev.off()
  
  
  # outlier detection, as suggedted by (Theis, 2019)
  
  outlier_detect <- function(a_vector, low_q=0.1, high_q=0.9){
    thresholds <- quantile(a_vector, c(low_q, high_q))
    (a_vector < thresholds[1]) | (a_vector > thresholds[2])
  }
  
  sce$outlier_counts <- sce$total_counts %>% outlier_detect()
  sce$outlier_features <- sce$total_features_by_counts %>% outlier_detect()
  sce$outlier_fcontrol <- sce$pct_counts_feature_control %>% outlier_detect()
  sce$outlier_mito <- sce$pct_counts_mitochondrial %>% outlier_detect()
  sce$outlier_hits <- sce$outlier_counts + sce$outlier_features +
    sce$outlier_fcontrol
  
  
  #saving per cell stats in a table
  per_cell_metric <- c("plate_position", "outlier_hits",
                        "outlier_counts" ,"outlier_features", "outlier_fcontrol",
                       "is_cell_control", "total_counts", "total_features_by_counts",
                       "pct_counts_mitochondrial", "pct_counts_ERCC", "pct_counts_in_top_100_features")
  
  per_cell_table <- as.data.frame(colData(sce)[, per_cell_metric])
  rownames(per_cell_table) <- rownames(per_cell_table) %>% 
    gsub("_[CATG]{8}-[CATG]{8}","",.)
  sample_names <- per_cell_table %>% rownames
  per_cell_table$sample_ID <- sample_names
  per_cell_table <- per_cell_table %>% select(sample_ID, everything())
  write.table(per_cell_table, paste0("Per_sample_key_metrics_",plate_id,".tsv"),
              sep='\t', row.names = F)
  
  triple_outliers <- which(per_cell_table$outlier_hits==3) %>% sample_names[.]
  double_outliers <- which(per_cell_table$outlier_hits==2) %>% sample_names[.]
  single_outliers <- which(per_cell_table$outlier_hits==1) %>% sample_names[.]
  
  save(plot_cumulative_dist, plot_ercc_scatter, plot_expr_vs_mean, plot_fctrl_scatter,
       plot_featurenum_vs_counts, plot_highest_expression, plot_mito_scatter,
       plot_plate_counts, plot_plate_mito, plot_plate_ercc,
       plot_pca, plot_tsne,
       per_cell_table, ercc_table, per_cell_metric,
       triple_outliers, double_outliers, single_outliers,
       plate_id, experiment_id,
       sce,
       file=paste0(plate_id,"qc_for_doc.Rdata"))
  
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
    outloc <- args[2]
    if(!file.exists(outloc)){
      warning(paste0("Output directory ('",outloc,"') does not exist!"),
              "\nCreating the directory...")
      system(paste("mkdir -p",outloc))
    }
  }
  if(length(args)>=3)
    plate_info <- args[3]
  if(length(args)>=4)
    sheetfile <- args[4]
  if(length(args)>=5)
    mtrds <- args[5]

  scqc_from_tsv(args[1], outloc, plate_info, sheetfile, mtrds)
}
