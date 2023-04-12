# based on /ei/cb/development/lany/scqc_reqs-0.2.dev/scqc_reqs-0.2.2/scqc_from_matrix.meta1.R
# current location: /ei/cb/development/lany/scqc_reqs-0.2.dev/scqc_reqs-0.2.1/scqc_from_matrix.meta1.v2.R
# updated: 
#    1. output QC_meanexp_vs_freq$plateID.pdf

library(dplyr)
library(scater)

#' Running single-cell QC on a tpm/estcount matrix
#'
#' @param tsv_location Expression matrix
#' @param output_location Where to store the QC files.
#' @param sheet_file Sample sheet with plate position, control status and sample ID information
#' If missing, plate position will try to be inferred from the .fastq file names
#' @param mt_file Mitochondrial genes .rds. If not provided, genes starting with "MT-" or
#' "mt-" will be considered mitochondrial (works only for human and mouse)
#'
#' @return None. Generates files (pdf plots, tsv tables).
#' @export
#'
#' @examples
# tsv_location="/ei/cb/analysis/CB-PPBFX-979_Iain_Macaulay_EI_CU_ENQ-5072_A_01/Analysis/run2.theleria_parva/scqc_reqs-0.1/quants_dir/est_countsCU5DAY0_matrix.tsv"


scqc_from_tsv <- function(tsv_location, pct_pseudo_txt, output_location="./",
                          sheet_file=NULL, mt_file=NULL) {
  plot_pct_TParva=F
  # test for GENANNO-552_Charlotte_Utting_EI_CU_ENQ-5476_A_01
  if (F){
    rundir='/Volumes/core-bioinformatics/analysis/CB-GENANNO-552_Charlotte_Utting_EI_CU_ENQ-5476_A_01/Analysis/scqc_reqs-0.2.4/run1'
    mt_file='/Volumes/core-bioinformatics/common/Databases/scqc/Bos_taurus/ARS-UCD1.2/MT_cow_vector.rds'
    
    # rundir='/ei/cb/analysis/CB-GENANNO-552_Charlotte_Utting_EI_CU_ENQ-5476_A_01/Analysis/scqc_reqs-0.2.4/run1'
    # mt_file='/ei/cb/common/Databases/scqc/Bos_taurus/ARS-UCD1.2/MT_cow_vector.rds'
    # tsv_location <- file.path(rundir,'quants_dir/est_countsCUB5709DAY69_matrix.tsv');
    tsv_location <- file.path(rundir,'quants_dir/all_plates.tsv');
    pct_pseudo_txt <- file.path(rundir, 'qc_dir/percent_pseudoaligned.txt');
    output_location <- file.path(rundir,'qc_dir');
    sheet_file <- file.path(rundir,'SampleSheetPIP-3109.v2.unix.csv');
    plot_pct_TParva=F
    is_test=F
  }
  # test for PPBFX979
  if (F) {
    tsv_location="/ei/cb/analysis/CB-PPBFX-979_Iain_Macaulay_EI_CU_ENQ-5072_A_01/Analysis/run2.theleria_parva/scqc_reqs-0.1/quants_dir/all_plates.tsv"
    pct_pseudo_txt <- "/ei/cb/analysis/CB-PPBFX-979_Iain_Macaulay_EI_CU_ENQ-5072_A_01/Analysis/run2.theleria_parva/scqc_reqs-0.1/qc_dir/percent_pseudoaligned.txt"
    output_location="/ei/cb/development/lany/CB-Iain_Macaulay_EI_CU_ENQ-5072_A_01/Analysis/scqc_reqs-0.2/test1/quants_dir/"
    sheet_file="/ei/cb/development/lany/CB-Iain_Macaulay_EI_CU_ENQ-5072_A_01/Analysis/scqc_reqs-0.2/test1/CowSeq_Datasheet.v3.unix.csv"
    mt_file="/ei/cb/common/Databases/scqc/Bos_taurus/ARS-UCD1.2/MT_cow_vector.rds"
    plot_pct_TParva=T
    is_test=T
    plate_id <- 'all'
  } else {
    plot_pct_TParva=F
    is_test=F
  }
  # test for PPBFX-996
  if (F) {
    tsv_location="/ei/cb/analysis/CB-PPBFX-996_Geoff_Mok_UOE_GM_ENQ-5130_A_01/Analysis/scqc_reqs-0.1/run1/quants_dir/all_plates.tsv"
    pct_pseudo_txt <- "/ei/cb/analysis/CB-PPBFX-996_Geoff_Mok_UOE_GM_ENQ-5130_A_01/Analysis/scqc_reqs-0.1/run1/qc_dir/percent_pseudoaligned.txt"
    output_location="/ei/cb/analysis/CB-PPBFX-996_Geoff_Mok_UOE_GM_ENQ-5130_A_01/Analysis/scqc_reqs-0.1/run1/quants_dir.test"
    sheet_file="/ei/cb/analysis/CB-PPBFX-996_Geoff_Mok_UOE_GM_ENQ-5130_A_01/Analysis/scqc_reqs-0.1/run1/SampleSheetPPBFX-996.v1.csv"
    mt_file="/ei/cb/common/Databases/scqc/Gallus_gallus/GRCg6a/MT_Gallus_gallus_vector.rds"
    plot_pct_TParva=F
    is_test=T
    plate_id <- 'all'
  } else {
    plot_pct_TParva=F
    is_test=F
  }
  # test for cow-seq scqc_req-0.2.1 debug, one plate
  if (F){
    tsv_location='/ei/cb/development/lany/CB-Iain_Macaulay_EI_CU_ENQ-5072_A_01/Analysis/scqc_reqs-0.2.1/test1/quants_dir/est_countsCU5DAY7_matrix.tsv'
    pct_pseudo_txt='/ei/cb/development/lany/CB-Iain_Macaulay_EI_CU_ENQ-5072_A_01/Analysis/scqc_reqs-0.2.1/test1/work/89/00f2c164b23995e467f193617be9ae/percent_pseudoaligned.txt'
    output_location='/ei/cb/development/lany/CB-Iain_Macaulay_EI_CU_ENQ-5072_A_01/Analysis/scqc_reqs-0.2.1/test1/qc_dir.1/'
    sheet_file='/ei/cb/development/lany/CB-Iain_Macaulay_EI_CU_ENQ-5072_A_01/Analysis/scqc_reqs-0.2/test1/CowSeq_Datasheet.v3.unix.csv'
    mt_file='/ei/cb/common/Databases/scqc/Bos_taurus/ARS-UCD1.2/MT_cow_vector.rds'
  }
  # test for cow-seq scqc_req-0.2.2 debug, all plates
  if (F){
    tsv_location='/ei/cb/development/lany/CB-Iain_Macaulay_EI_CU_ENQ-5072_A_01/Analysis/scqc_reqs-0.2.1/test1/quants_dir/all_plates.tsv'
    pct_pseudo_txt='/ei/cb/development/lany/CB-Iain_Macaulay_EI_CU_ENQ-5072_A_01/Analysis/scqc_reqs-0.2.1/test1/qc_dir/percent_pseudoaligned.txt'
    output_location='/ei/cb/development/lany/CB-Iain_Macaulay_EI_CU_ENQ-5072_A_01/Analysis/scqc_reqs-0.2.1/test1/qc_dir.2/'
    sheet_file='/ei/cb/development/lany/CB-Iain_Macaulay_EI_CU_ENQ-5072_A_01/Analysis/scqc_reqs-0.2/test1/CowSeq_Datasheet.v3.unix.csv'
    mt_file='/ei/cb/common/Databases/scqc/Bos_taurus/ARS-UCD1.2/MT_cow_vector.rds'
    plot_pct_TParva=F
    is_test=F
  }
  # test for cowseq2 hpc
  if (F){
    tsv_location='/ei/cb/analysis/CB-GENANNO-515_Charlotte_Utting_EI_CU_ENQ-5187_A_01/Analysis/scqc_reqs-0.2.1/run1/quants_dir/all_plates.tsv'
    pct_pseudo_txt='/ei/cb/analysis/CB-GENANNO-515_Charlotte_Utting_EI_CU_ENQ-5187_A_01/Analysis/scqc_reqs-0.2.1/run1/qc_dir/percent_pseudoaligned.txt'
    output_location='/ei/cb/analysis/CB-GENANNO-515_Charlotte_Utting_EI_CU_ENQ-5187_A_01/Analysis/scqc_reqs-0.2.1/run1/qc_dir/'
    sheet_file='/ei/cb/analysis/CB-GENANNO-515_Charlotte_Utting_EI_CU_ENQ-5187_A_01/Analysis/scqc_reqs-0.2.1/run1/cowseq_Batch2-data.PSEQ-2345.2.unix.csv'
    mt_file='/ei/cb/common/Databases/scqc/Bos_taurus/ARS-UCD1.2/MT_cow_vector.rds'
  }
  # test for cowseq2 local 1 plate
  if (F){
    rundir='/Volumes/core-bioinformatics/analysis/CB-GENANNO-515_Charlotte_Utting_EI_CU_ENQ-5187_A_01/Analysis/scqc_reqs-0.2.2/run1/'
    tsv_location <- file.path(rundir,'quants_dir/est_countsCUB25DAY0_matrix.tsv');
    pct_pseudo_txt <- file.path(rundir, 'qc_dir/percent_pseudoaligned.txt');
    output_location <- file.path(rundir,'qc_dir/');
    sheet_file <- file.path(rundir,'cowseq_Batch2-data.PSEQ-2345.2.unix.csv');
    mt_file='/Volumes/core-bioinformatics/common/Databases/scqc/Bos_taurus/ARS-UCD1.2/MT_cow_vector.rds'
    plot_pct_TParva=T
    is_test=F
  }
  # test for cowseq2 local all plates
  if (F){
    rundir='/Volumes/core-bioinformatics/analysis/CB-GENANNO-515_Charlotte_Utting_EI_CU_ENQ-5187_A_01/Analysis/scqc_reqs-0.2.2/run1/'
    tsv_location <- file.path(rundir,'quants_dir/all_plates.tsv');
    pct_pseudo_txt <- file.path(rundir, 'qc_dir/percent_pseudoaligned.txt');
    output_location <- file.path(rundir,'qc_dir/qc_all');
    sheet_file <- file.path(rundir,'cowseq_Batch2-data.PSEQ-2345.2.unix.csv');
    mt_file='/Volumes/core-bioinformatics/common/Databases/scqc/Bos_taurus/ARS-UCD1.2/MT_cow_vector.rds'
    plot_pct_TParva=T
    is_test=F
  }
  # test for chicken2 1 plate
  if (F){
    rundir='/Volumes/core-bioinformatics/analysis/CB-GENANNO-516_Gi_Fay_Mok_UOE_GM_ENQ-5204_A_01/Analysis/scqc_reqs-0.2.2/run1/'
    tsv_location <- file.path(rundir,'quants_dir/est_countsGMPlate16w_matrix.tsv');
    pct_pseudo_txt <- file.path(rundir, 'qc_dir/percent_pseudoaligned.txt');
    output_location <- file.path(rundir,'qc_dir');
    sheet_file <- file.path(rundir,'AS_SampleSheetPSEQ-2347.2.unix.csv');
    mt_file='/Volumes/core-bioinformatics/common/Databases/scqc/Gallus_gallus/GRCg6a/MT_Gallus_gallus_vector.rds'
    plot_pct_TParva=F
    is_test=F
  }
  # test for human control 1 plate, no meta data
  if (F){
    rundir='/ei/cb/development/lany/CB-GENANNO-523_Karim_Gharbi_EI_KG_ENQ-3708/Analysis/scqc_reqs-0.2.2/run1'
    tsv_location <- file.path(rundir,'quants_dir/est_countsR0121-P0002_matrix.tsv');
    pct_pseudo_txt <- file.path(rundir, 'qc_dir/percent_pseudoaligned.txt');
    output_location <- file.path(rundir,'qc_dir');
    sheet_file <- file.path(rundir,'SampleSheetPSEQ-2015.v2.unix.csv');
    mt_file='/ei/cb/common/Databases/scqc/Homo_sapiens.with.ERCC/GRCh38_and_ERCC92/MT_GRCh38_vector.rds'
    plot_pct_TParva=F
    is_test=F
  }
  # test for human control 1 plate, with meta data
  if (F){
    # rundir='/ei/cb/development/lany/CB-GENANNO-523_Karim_Gharbi_EI_KG_ENQ-3708/Analysis/scqc_reqs-0.2.2/run2.metadata'
    rundir="/Volumes/core-bioinformatics/development/lany/CB-GENANNO-523_Karim_Gharbi_EI_KG_ENQ-3708/Analysis/scqc_reqs-0.2.2/run2.metadata"
    tsv_location <- file.path(rundir,'quants_dir/est_countsR0121-P0001_matrix.tsv');
    pct_pseudo_txt <- file.path(rundir, 'qc_dir/percent_pseudoaligned.txt');
    output_location <- file.path(rundir,'qc_dir');
    sheet_file <- file.path(rundir,'SampleSheetPSEQ-2015.metadata.v3.unix.csv');
    # mt_file='/ei/cb/common/Databases/scqc/Homo_sapiens.with.ERCC/GRCh38_and_ERCC92/MT_GRCh38_vector.rds'
    mt_file='/Volumes/core-bioinformatics/common/Databases/scqc/Homo_sapiens.with.ERCC/GRCh38_and_ERCC92/MT_GRCh38_vector.rds'
    plot_pct_TParva=F
    is_test=F
  }
  # test for PBFX-1067_Claudia_Ribeiro_de_Almeida, meta1
  if (F){
    rundir='/ei/cb/analysis/CB-PPBFX-1067_Claudia_Ribeiro_de_Almeida_BI_CR_ENQ-5220_A_01_Smart-Seq2_QC/Analysis/scqc_reqs-0.2.2/run2'
    tsv_location <- file.path(rundir,'quants_dir/est_countsR0881-P0004_matrix.tsv');
    # tsv_location <- file.path(rundir,'quants_dir/all_plates.tsv');
    pct_pseudo_txt <- file.path(rundir, 'qc_dir/percent_pseudoaligned.txt');
    output_location <- file.path(rundir,'qc_dir');
    sheet_file <- file.path(rundir,'SampleSheetPSEQ-2422_v3_SingleCellQC.2.unix.csv');
    mt_file='/ei/cb/common/Databases/scqc/Mus_musculus/GRCm39/MT.Mus_musculus.GRCm39.rds'
    plot_pct_TParva=F
    is_test=F
  }
  # test for GENANNO-543_Stuart_Rushworth_UOE_SR_ENQ-5327, meta1
  if (F){
    rundir='/ei/cb/analysis/CB-GENANNO-543_Stuart_Rushworth_UOE_SR_ENQ-5327_A_01/Analysis/scqc_reqs-0.2.2/run2'
    #tsv_location <- file.path(rundir,'quants_dir/est_countsR0881-P0004_matrix.tsv');
    # tsv_location <- file.path(rundir,'quants_dir/all_plates.tsv');
    pct_pseudo_txt <- file.path(rundir, 'qc_dir/percent_pseudoaligned.txt');
    output_location <- file.path(rundir,'qc_dir');
    sheet_file <- file.path(rundir,'META_DATA_CU_EW_v3.unix.csv');
    mt_file='/ei/cb/common/Databases/scqc/Mus_musculus/GRCm39/MT.Mus_musculus.GRCm39.rds'
    plot_pct_TParva=F
    is_test=F
  }
  # test for png output on GENANNO-543_Stuart_Rushworth_UOE_SR_ENQ-5327, meta1
  if (F){
    rundir='/ei/cb/analysis/CB-GENANNO-543_Stuart_Rushworth_UOE_SR_ENQ-5327_A_01/Analysis/scqc_reqs-0.2.2/run1/'
    tsv_location <- file.path(rundir,'quants_dir/est_countsLTHSCO1_matrix.tsv');
    # tsv_location <- file.path(rundir,'quants_dir/all_plates.tsv');
    pct_pseudo_txt <- file.path(rundir, 'qc_dir/percent_pseudoaligned.txt');
    output_location <- file.path(rundir,'qc_dir');
    sheet_file <- file.path(rundir,'META_DATA_CU_EW_v3.unix.csv');
    mt_file='/ei/cb/common/Databases/scqc/Mus_musculus/GRCm39/MT.Mus_musculus.GRCm39.rds'
    plot_pct_TParva=F
    is_test=F
  }
  # test for Charlotte_Utting_MiSeq, meta1
  if (F){
    rundir='/ei/cb/development/lany/Charlotte_Utting_MiSeq/Analysis/scqc_reqs-0.2.2/run2'
    tsv_location <- file.path(rundir,'quants_dir/est_countsP0001_matrix.tsv');
    # tsv_location <- file.path(rundir,'quants_dir/all_plates.tsv');
    pct_pseudo_txt <- file.path(rundir, 'qc_dir/percent_pseudoaligned.txt');
    output_location <- file.path(rundir,'qc_dir');
    sheet_file <- file.path(rundir,'SampleSheetUsed.v2.unix.csv');
    mt_file='/ei/cb/common/Databases/scqc/Mus_musculus/GRCm39/MT.Mus_musculus.GRCm39.rds'
    plot_pct_TParva=F
    is_test=F
  }
  # ===== 1. read expression matrix =====
  #tsv will be est_counts or tpm matrix transcript X sample
  #but counts still used as label in the SCE object
  counts_df <- read.table(tsv_location, header = T, check.names = F)
  fastqname <- colnames(counts_df)
  # experiment_id <- gsub("/[A-z|_|\\.]+$","",tsv_location) %>% 
  #   gsub("^.*/","",.)
  plate_id <- "1"
  
  # 1.1: using the sheet_file if provided
  if(!is.null(sheet_file)){
    samplesheet <- read.csv(sheet_file)
    # ==== check what columns are missing from sample sheet ====
    column_to_check=c('Sample_Plate', 'Sample_Well', 'control', 'number_of_cells','meta_1','meta_2','ercc')
    samplesheet.check <- data.frame(t(is.element(column_to_check,colnames(samplesheet))))
    colnames(samplesheet.check) <- column_to_check
    # > samplesheet.check[1,]
    #   Sample_Plate Sample_Well control number_of_cells meta_1 meta_2  ercc
    # 1         TRUE        TRUE   FALSE           FALSE  FALSE  FALSE FALSE
    
    # check if above columns are fully entried (as T)
    # colSums(is.na.data.frame(samplesheet))==0
    # >     colSums(is.na.data.frame(samplesheet))==0
    # Lane                Sample_ID             Sample_Name            Sample_Plate 
    # TRUE                    TRUE                    TRUE                    TRUE 
    # Sample_Well             i7_index_ID            Index             i5_index_ID 
    # TRUE                   FALSE                    TRUE                   FALSE 
    # Index2          Sample_Project             LibraryType                 TaxonID 
    # TRUE                    TRUE                    TRUE                    TRUE 
    # ScientificName DataQCPipeline                   UUID         DataAccessGroup 
    # TRUE                   FALSE                   FALSE                    TRUE 
    # unique_sample_id_suffix 
    # TRUE 
    # 
    # make sure a column has all entries as na
    
    samplesheet.check[1,] <- samplesheet.check[1,]&
      (colSums(is.na.data.frame(samplesheet))<nrow(samplesheet))[colnames(samplesheet.check)]
    # > samplesheet.check
    #   Sample_Plate Sample_Well control number_of_cells meta_1 meta_2  ercc
    # 1         TRUE        TRUE   FALSE           FALSE  FALSE  FALSE FALSE
    if (!samplesheet.check$Sample_Well) stop('SampleSheet ERROR: "Sample_Well" column missing/incomplete!')
    if (!samplesheet.check$Sample_Plate) stop('SampleSheet ERROR: "Sample_Plate" column missing/incomplete!')
    
    # internally derive unique_sample_id_suffix. use it to link fastqs
    samplesheet$unique_sample_id_suffix <- paste0(samplesheet$Sample_Name,'_',samplesheet$Sample_ID)
    row_ids <- sapply(samplesheet$unique_sample_id_suffix,
           function(id) grepl(paste0(id, "_"), colnames(counts_df))) %>% apply(1,which)
    # > samplesheet[row_ids[1:3],]
    #     Lane Sample_ID             Sample_Name Sample_Plate Sample_Well i7_index_ID    Index
    # 281    1 A06956_1 R0121-S0097_2V2aA01s481  R0121-P0002          A1          NA TAAGGCGA
    # 282    1 A06957_1 R0121-S0098_2V2aB01s482  R0121-P0002          B1          NA TAAGGCGA
    # 283    1 A06958_1 R0121-S0099_2V2aC01s483  R0121-P0002          C1          NA TAAGGCGA
    #     i5_index_ID   Index2 Sample_Project            LibraryType TaxonID ScientificName
    # 281          NA CTAGTCGA       PIP-2406 Illuminapremadelibrary    9606   Homo sapiens
    # 282          NA AGCTAGAA       PIP-2406 Illuminapremadelibrary    9606   Homo sapiens
    # 283          NA ACTCTAGG       PIP-2406 Illuminapremadelibrary    9606   Homo sapiens
    #     DataQCPipeline UUID DataAccessGroup          unique_sample_id_suffix
    # 281             NA   NA        Internal R0121-S0097_2V2aA01s481_A06956_1
    # 282             NA   NA        Internal R0121-S0098_2V2aB01s482_A06957_1
    # 283             NA   NA        Internal R0121-S0099_2V2aC01s483_A06958_1

    # experiment_id <- samplesheet[row_ids,]$experiment %>% unique
    colnames(counts_df) <- samplesheet[row_ids,]$unique_sample_id_suffix
    # > counts_df[1:2,1:4]
    #                            R0121-S0097_2V2aA01s481_A06956_1 R0121-S0098_2V2aB01s482_A06957_1
    # transcript:ENST00000456328                                0                                0
    # transcript:ENST00000450305                                0                                0
    #                            R0121-S0099_2V2aC01s483_A06958_1 R0121-S0100_2V2aD01s484_A06959_1
    # transcript:ENST00000456328                                0                                0
    # transcript:ENST00000450305
    plate_id <- as.character(unique(samplesheet[row_ids,"Sample_Plate"]))
  }
  
  if(length(plate_id)>1){
    plate_id <- "all"
  }
  
  if (!samplesheet.check$control) {
    inds_cellctrl <- numeric(0)
  } else {
    inds_cellctrl <- which(samplesheet[row_ids,]$control==T)
  }
  #feature and cell controls
  inds_ercc <- counts_df %>% rownames() %>% grep("^ERCC",.)
  inds_mito <- counts_df %>% rownames() %>% grep("\\|MT-",., ignore.case = T) #only for mouse?
  
  
  long_names <- counts_df[-inds_ercc,] %>% rownames()
  #takes a while, for transcripts with form "geneID|XXXX", this will remove any text after "|"
  short_names <- sapply(long_names, function(s) strsplit(s,"\\|") %>% .[[1]] %>% .[1]) %>% as.vector 
  rownames(counts_df)[-inds_ercc] <- short_names

  mito_names <- long_names[inds_mito] %>% strsplit("\\|") %>% sapply(function(x) x[5]) %>% as.vector
  rownames(counts_df)[inds_mito] <- mito_names
  
  if(!is.null(mt_file)){
    mito_names <- readRDS(mt_file)
    inds_mito <- counts_df %>% rownames %>%
      match(mito_names, .) %>% na.omit()
  }
  
  if (is_test){
    save(counts_df, inds_ercc, inds_mito, inds_cellctrl, row_ids, cell_number, well_id,plate_id,plate_id.allcells,
         experiment_id,samplesheet, file=file.path(output_location,sprintf('%s_countsdf.Rdata',plate_id)))
  }
  
  # ==== single cell analysis starts here =====
  if (is_test){
    load(file.path(output_location,sprintf('%s_countsdf.Rdata','all')), verbose = T)
  }
  sce <- SingleCellExperiment(assays=list(counts=counts_df %>% as.matrix))
  sce$unique_sample_id_suffix <- colnames(counts_df)
  sce$fastqname <- fastqname
  # sce$Cell <- colnames(counts(sce)) # this is not used anywere
  sce$control <- NULL
  if (samplesheet.check$control) sce$control<-samplesheet[row_ids,]$control
  sce$plate_position <- samplesheet[row_ids,]$Sample_Well

  sce$meta1 <- NULL
  if (samplesheet.check$meta_1) sce$meta1<-as.factor(samplesheet[row_ids,]$meta_1)
  sce$meta2 <- NULL
  if (samplesheet.check$meta_2) sce$meta2<-as.factor(samplesheet[row_ids,]$meta_2)
  # for case number_of_cells is 0,1
  # if (sum(unique(samplesheet[row_ids,]$number_of_cells)==c(0,1))==2){
  #   sce$nb_cells<-factor(samplesheet[row_ids,]$number_of_cells, 
  #                        levels = c('1','0'))
  # } else {
  sce$nb_cells <- 1
  if (samplesheet.check$number_of_cells) sce$nb_cells<-as.factor(samplesheet[row_ids,]$number_of_cells)
  # }
  
  sce$plate_id <- as.factor(samplesheet[row_ids,]$Sample_Plate)
  # column meta1_numeric and meta2_numeric: currently not implemented
  if (F){
    sce$meta1_numeric <- NULL
    if (!is.null(samplesheet$meta1_numeric)) sce$meta1_numeric<-samplesheet[row_ids,]$meta_1_numeric
    sce$meta2_numeric <- NULL
    if (!is.null(samplesheet$meta2_numeric)) sce$meta2_numeric<-samplesheet[row_ids,]$meta_2_numeric
    
  }

  # sce <- calculateQCMetrics(sce,
  #                           feature_controls = list(ERCC = inds_ercc, mitochondrial = inds_mito),
  #                           cell_controls = list(Controls = inds_cellctrl) 
  #                             )
  sce <- addPerCellQC(sce, 
                               subsets=list(ERCC = inds_ercc, mitochondrial = inds_mito))  
  sce <- addPerFeatureQC(sce, subsets=list(cellctrl = inds_cellctrl))
  rowData(sce)$total_exp_cells <-  scater::nexprs(sce, byrow=T, detection_limit=0)
  rowData(sce)$pct_exp_cells <- rowData(sce)$total_exp_cells / ncol(sce) *100
  rowData(sce)$is_MT <- F
  if (length(inds_mito)>0) {
    rowData(sce)$is_MT[inds_mito] <- T
  }
  rowData(sce)$is_ERCC <- F
  if (length(inds_ercc)>0) {
    rowData(sce)$is_ERCC[inds_ercc] <- T
  }
  
  # > head(colData(sce))
  # DataFrame with 6 rows and 18 columns
  #                                           unique_sample_id_suffix   plate_position  nb_cells    plate_id              sum  detected ...
  # R0121-S0097_2V2aA01s481_A06956_1 R0121-S0097_2V2aA01s481_A06956_1               A1         1 R0121-P0002 256.999995494704       247
  # R0121-S0098_2V2aB01s482_A06957_1 R0121-S0098_2V2aB01s482_A06957_1               B1         1 R0121-P0002 17891.9997868149      7229
  # R0121-S0099_2V2aC01s483_A06958_1 R0121-S0099_2V2aC01s483_A06958_1               C1         1 R0121-P0002 14966.9988926191      6755
  # R0121-S0100_2V2aD01s484_A06959_1 R0121-S0100_2V2aD01s484_A06959_1               D1         1 R0121-P0002 22724.9999667178      8832
  # R0121-S0101_2V2aE01s485_A06960_1 R0121-S0101_2V2aE01s485_A06960_1               E1         1 R0121-P0002 22625.0000854361      8159
  # R0121-S0102_2V2aF01s486_A06961_1 R0121-S0102_2V2aF01s486_A06961_1               F1         1 R0121-P0002 30711.0014675777     10666

  # specific to B.Taurus with T.Parva
  if (plot_pct_TParva){
    # specific to B.Taurus with T.Parva
    rowData(sce)$species <- rep(NA,nrow(rowData(sce)))
    rowData(sce)$species[grep('TpMuguga',rownames(rowData(sce)))] <- 'T_Parva'
    rowData(sce)$species[grep('transcript:ENSBTAT',rownames(rowData(sce)))] <- 'B_Taurus'
    # calculate per cell %count from each of TParva and BTaurus
    colData(sce)$pct_counts_in_TParva <- rep(0.0000, nrow(colData(sce)))
    colData(sce)$total_counts_in_TParva <- rep(0, nrow(colData(sce)))
    
    species_row <- rowData(sce)$species=='T_Parva'
    pct_count_in_species <- function(one_column, species_row){
      return(sum(one_column[species_row])/sum(one_column))
    }
    colData(sce)$pct_counts_in_TParva <- apply(X=counts(sce), MARGIN = 2, 
                                               FUN = pct_count_in_species, species_row=species_row)
    total_count_in_species <- function(one_column, species_row){
      return(sum(one_column[species_row]))
    }
    colData(sce)$total_counts_in_TParva <- apply(X=counts(sce), MARGIN = 2, 
                                               FUN = total_count_in_species, species_row=species_row)
  }
  # rm(counts_df)
  setwd(output_location)
  
  # ===== 1. QC_scatterplots_X.pdf ====
  cellctrl_flag <- NULL
  if(length(inds_cellctrl)>0){
    cellctrl_flag <- "control"
  }
  nbcell_flag <- NULL
  if(length(unique(sce$nb_cells))>1){
    nbcell_flag <- "nb_cells"
  }
  # for case: nbcells=(0,1) with 0 as cell control, we want 1 to have small marker size.
  cellsize_flip_flag=F
  nbcell_markersize=NULL
  if (F){
    # the case when number_of_cells=(0,1)
    if (length(unique(samplesheet[row_ids,]$number_of_cells))==2){
      if (sum(sort(unique(samplesheet[row_ids,]$number_of_cells))==c(0,1))==2) cellsize_flip_flag=T
    }
  }
  if (F){
    if (length(unique(samplesheet[row_ids,]$number_of_cells))>1){
      # a case: assign the most frequent number to smallest size
      t <- samplesheet[row_ids,'number_of_cells'] %>% table %>% sort(.,decreasing = T) %>% as.data.frame 
      colnames(t) <- c('number_of_cells','freq')
      t$number_of_cells <- as.numeric(as.character(t$number_of_cells))
      # > t
      # number_of_cells freq
      # 1               1   82
      # 2               0    6
      # 3               2    6
      # 4              20    1
      # 5              50    1
      cellsize_flip_flag=T
      t$markersize <- seq(2,2*nrow(t),by=2)
      t <- t[order(t$number_of_cells),]
      nbcell_markersize <- t$markersize
    }
  }
  
  if (length(unique(samplesheet[row_ids,]$number_of_cells))>2){
    if (is.element(1, unique(samplesheet[row_ids,]$number_of_cells))){
      cellsize_flip_flag=T
      t <- data.frame(number_of_cells=sort(unique(samplesheet[row_ids,]$number_of_cells)),
                      markersize=0)
      markersize <- seq(2,2*nrow(t),by=2)
      t$markersize[t$number_of_cells==1]=2
      t$markersize[t$number_of_cells!=1]=markersize[markersize!=2]
      nbcell_markersize <- t$markersize
    }
  }
  
  meta1_flag <- NULL
  if (length(unique(sce$meta1))>1) meta1_flag <- 'meta1'
  meta2_flag <- NULL
  if (length(unique(sce$meta2))>1) meta2_flag <- 'meta2'
  
  MT_flag <- NULL
  if(length(inds_mito)>0) MT_flag <- "is_MT"
  
  ERCC_flag <- NULL
  if(length(inds_ercc)>0) ERCC_flag <- "is_ERCC"
  
  pdf(paste0("QC_scatterplots_",plate_id,".pdf"), pointsize = 12)
  
  # plot1: log_count vs total feature
  # optimal display:
  call_plotColData <- function(sce=NULL, x = "sum",y = "detected",
                               colour_by = NULL, size_by= NULL, shape_by=NULL, cellsize_flip_flag=F,
                               nbcell_markersize=NULL
                               ) {
    
    p <- plotColData(sce, x = x,y = y, colour_by = colour_by, size_by= size_by, shape_by=shape_by)
    p <- p+scale_x_log10()
    if (cellsize_flip_flag) p <- p + scale_size_manual(values=nbcell_markersize)
    p <- p+theme(legend.position = "right") + xlab("log_total_counts") + ylab("total_features") +
      labs(title="Total features ~ total counts (log scale)")
    return(p)
  }
  plot_featurenum_vs_counts <- list()
  pidx <- 1
  size_by= nbcell_flag
  if (!is.null(cellctrl_flag)){
    colour_by = cellctrl_flag
    shape_by=NULL
    if (!is.null(meta1_flag)){
      shape_by=meta1_flag
      plot_featurenum_vs_counts[[pidx]] <- 
        call_plotColData(sce, colour_by = colour_by, size_by= size_by,cellsize_flip_flag=cellsize_flip_flag, 
                                nbcell_markersize = nbcell_markersize,shape_by=meta1_flag)
      
      if (F){# this is to test the plot, make sure all text fit inside figure
        p <- plotColData(sce, x = "sum",y = "detected",
                         colour_by = colour_by, size_by= size_by, shape_by=shape_by)
        p <- p+scale_x_log10()
        if (cellsize_flip_flag) p <- p + scale_size_manual(values=c(4,2,6,8,10))
        p <- p+theme(legend.position = "right") + xlab("log_total_counts") + ylab("total_features") +
          labs(title="Total features ~ total counts (log scale)")
      }
      plot(plot_featurenum_vs_counts[[pidx]])
      pidx <- pidx+1
    }
    if (!is.null(meta2_flag)){
      shape_by=meta2_flag
      plot_featurenum_vs_counts[[pidx]] <-call_plotColData(sce, colour_by = colour_by, size_by= size_by,
                               cellsize_flip_flag=cellsize_flip_flag, nbcell_markersize = nbcell_markersize,shape_by=meta2_flag)
      plot(plot_featurenum_vs_counts[[pidx]])
      pidx <- pidx+1
    }
    if(is.null(meta1_flag) &  is.null(meta2_flag)){
      plot_featurenum_vs_counts[[pidx]] <-call_plotColData(sce, colour_by = colour_by, size_by= size_by,
                                                           cellsize_flip_flag=cellsize_flip_flag, nbcell_markersize = nbcell_markersize,shape_by=NULL)
      plot(plot_featurenum_vs_counts[[pidx]])
      pidx <- pidx+1
    }
  } else {
    shape_by=NULL
    colour_by=NULL
    if (!is.null(meta1_flag)){
      plot_featurenum_vs_counts[[pidx]] <-call_plotColData(sce, colour_by=meta1_flag, size_by=size_by, shape_by=shape_by)
      plot(plot_featurenum_vs_counts[[pidx]])
      pidx <- pidx+1}
    if (!is.null(meta2_flag)){
      plot_featurenum_vs_counts[[pidx]] <-call_plotColData(sce, colour_by=meta2_flag, size_by=size_by, shape_by=shape_by)
      plot(plot_featurenum_vs_counts[[pidx]])
      pidx <- pidx+1}
    if(is.null(meta1_flag) &  is.null(meta2_flag)){
      plot_featurenum_vs_counts[[pidx]] <-call_plotColData(sce, colour_by=NULL, size_by=size_by, shape_by=shape_by)
      plot(plot_featurenum_vs_counts[[pidx]])
      pidx <- pidx+1}
  }
  
  if(plot_pct_TParva & !is.null(meta1_flag)){
    plot_featurenum_vs_counts[[pidx]] <- call_plotColData(sce, x='pct_counts_in_TParva', y='detected', 
                            colour_by=meta1_flag, size_by= nbcell_flag,cellsize_flip_flag=cellsize_flip_flag, shape_by=meta1_flag)+
      theme(legend.position = "top") + xlab("pct_counts_in_TParva (0,1)") + ylab("total_features") +
        labs(title="Total features ~ pct counts (log scale)")
    plot(plot_featurenum_vs_counts[[pidx]])
    pidx <- pidx+1

    plot_featurenum_vs_counts[[pidx]]<- call_plotColData(sce, x='total_counts_in_TParva', y='detected', 
                            colour_by=meta1_flag, size_by= nbcell_flag,cellsize_flip_flag=cellsize_flip_flag, shape_by=meta1_flag) + 
      theme(legend.position = "top") + xlab("total_counts_in_TParva (0,1)") + ylab("total_features") +
        labs(title="Total features ~ total counts TParva (log scale)")
    plot(plot_featurenum_vs_counts[[pidx]])
    pidx <- pidx+1
  }
  
  # plot2.
  plot_cumulative_dist <- plotScater(sce, nfeatures = 300, exprs_values = "counts",
                                     colour_by = cellctrl_flag) 
  plot(plot_cumulative_dist)
  
  if (F){
    idx_row=288 # 
    # sce$de
    
  }
  # plot3: mean expression vs frequency expressed, as part of the scatter pdf file
  if (F){
    plot_expr_vs_mean <- plotRowData(sce, x='mean', y='pct_exp_cells', by_exprs_values='counts', colour_by=ERCC_flag)+scale_x_log10()+
      ylab('percentage of expressing cells')+xlab('mean expression')+
      geom_hline(yintercept = 50, linetype='dashed')+geom_vline(xintercept = mean(assay(sce,'counts')),linetype='dashed')+
      geom_text(label=sprintf('%d genes are expressed in at least 50%% of cells', sum(rowData(sce)$pct_exp_cells>=50)), 
                y=40, x=10^(-05), size=4)+
      geom_text(label=sprintf('%d genes are expressed in at least 25%% of cells', sum(rowData(sce)$pct_exp_cells>=25)), 
                y=20, x=10^(-05), size=4)+
      geom_text(label=sprintf('%d genes are expressed in at least 1 cell', sum(rowData(sce)$detected>0)), 
                y=60, x=10^(-05), size=4)
    plot(plot_expr_vs_mean)
  }
  # plot4: total count across genes
  plot_highest_expression <- plotHighestExprs(sce, exprs_values = "counts", colour_cells_by=cellctrl_flag,
                                              as_percentage=T)+xlab('% of counts')
  plot(plot_highest_expression)
 
  #% plot 5 and the rest: any feature control
  if (!is.null(cellctrl_flag)){
    # MT
    plot_mito_scatter <- plotColData(sce, x = "detected", y = "subsets_mitochondrial_percent",
                                     colour_by = cellctrl_flag, size_by=nbcell_flag) + theme(legend.position = "right") + xlab("total_features") +
      labs(title=sprintf("Out of %d transcripts, Mitochondrial expression percentage ~ total features", nrow(counts(sce))))
    if (cellsize_flip_flag) plot_mito_scatter <- plot_mito_scatter + scale_size_manual(values=nbcell_markersize)
      
    plot(plot_mito_scatter)
    # ERCC
    plot_ercc_scatter <- plotColData(sce, x = "detected", y = "subsets_ERCC_percent",
                                     colour_by = cellctrl_flag, size_by=nbcell_flag) +  theme(legend.position = "right") + xlab("total_features") +
      labs(title="Spike-in expression percentage ~ total features")
    if (cellsize_flip_flag) plot_ercc_scatter <- plot_ercc_scatter + scale_size_manual(values=nbcell_markersize)
    plot(plot_ercc_scatter)
  } else {
    # MT
    plot_mito_scatter <- plotColData(sce, x = "detected", y = "subsets_mitochondrial_percent",
                                     colour_by = nbcell_flag, size_by=nbcell_flag) + theme(legend.position = "right") + xlab("total_features") +
      labs(title=sprintf("Out of %d transcripts, Mitochondrial expression percentage ~ total features", nrow(counts(sce))))
    if (cellsize_flip_flag) plot_mito_scatter <- plot_mito_scatter + scale_size_manual(values=nbcell_markersize)
    plot(plot_mito_scatter)
    # ERCC
    plot_ercc_scatter <- plotColData(sce, x = "detected", y = "subsets_ERCC_percent",
                                     colour_by = nbcell_flag, size_by=nbcell_flag) +  theme(legend.position = "right") + xlab("total_features") +
      labs(title="Spike-in expression percentage ~ total features")
    if (cellsize_flip_flag) plot_ercc_scatter <- plot_ercc_scatter + scale_size_manual(values=nbcell_markersize)
    plot(plot_ercc_scatter)
  }
  dev.off()
  
  # plot3: mean expression vs frequency expressed. save as pdf - will be converted to png using ghostscript.
  if (F){
    pdf(paste0("/ei/cb/analysis/CB-GENANNO-543_Stuart_Rushworth_UOE_SR_ENQ-5327_A_01/Analysis/scqc_reqs-0.2.2/run1/qc_dir.test/",
      "QC_meanexp_vs_freq",plate_id,".pdf"), pointsize = 8)
    plot_expr_vs_mean <- plotRowData(sce, x='mean', y='pct_exp_cells', by_exprs_values='counts', colour_by=ERCC_flag)+scale_x_log10()+
      ylab('percentage of expressing cells')+xlab('mean expression')+
      geom_hline(yintercept = 50, linetype='dashed')+geom_vline(xintercept = mean(assay(sce,'counts')),linetype='dashed')+
      geom_text(label=sprintf('%d genes are expressed in at least 50%% of cells', sum(rowData(sce)$pct_exp_cells>=50)), 
                y=40, x=10^(-05), size=4)+
      geom_text(label=sprintf('%d genes are expressed in at least 25%% of cells', sum(rowData(sce)$pct_exp_cells>=25)), 
                y=20, x=10^(-05), size=4)+
      geom_text(label=sprintf('%d genes are expressed in at least 1 cell', sum(rowData(sce)$detected>0)), 
                y=60, x=10^(-05), size=4)
    plot(plot_expr_vs_mean)
    dev.off()
    
  }
  # plot 3: as above, now output single pdf file - this will be converted to png via ghostscript
  if (T){
    pdf(paste0("QC_meanexp_vs_freq",plate_id,".pdf"), pointsize = 8)
    
    plot_expr_vs_mean <- plotRowData(sce, x='mean', y='pct_exp_cells', by_exprs_values='counts', colour_by=ERCC_flag)+scale_x_log10()+
      ylab('percentage of expressing cells')+xlab('mean expression')+
      geom_hline(yintercept = 50, linetype='dashed')+geom_vline(xintercept = mean(assay(sce,'counts')),linetype='dashed')+
      geom_text(label=sprintf('%d genes are expressed in at least 50%% of cells', sum(rowData(sce)$pct_exp_cells>=50)), 
                y=40, x=10^(-05), size=4)+
      geom_text(label=sprintf('%d genes are expressed in at least 25%% of cells', sum(rowData(sce)$pct_exp_cells>=25)), 
                y=20, x=10^(-05), size=4)+
      geom_text(label=sprintf('%d genes are expressed in at least 1 cell', sum(rowData(sce)$detected>0)), 
                y=60, x=10^(-05), size=4)
    plot(plot_expr_vs_mean)
    dev.off()
    
  }
  
  # plot3: mean expression vs frequency expressed. save using png function: THIS ONLY WORK ON RSTUDIO WITH GRAPHIC SUPPORT.
  if (F){
    #options(bitmapType='Xlib') doesn't work
    #options(bitmapType='cairo')doesn't work
    #options(bitmapType = 'cairo', device = 'png')doesn't work
    if (F){
      # test png output function. save to qc_dir.test
      png(paste0('/ei/cb/analysis/CB-GENANNO-543_Stuart_Rushworth_UOE_SR_ENQ-5327_A_01/Analysis/scqc_reqs-0.2.2/run1/qc_dir.test/',
                 "QC_meanexp_vs_freq",plate_id,".png"), pointsize = 12, width = 720, height = 720)
      
    }
    png(paste0("QC_meanexp_vs_freq",plate_id,".png"), pointsize = 12, width = 720, height = 720)
    # plot_expr_vs_mean <-  plotExprsFreqVsMean(sce, exprs_values = "counts", show_smooth = F)
    # function plotExprsFreqVsMean is discontinuous, so make our own
    plot_expr_vs_mean <- plotRowData(sce, x='mean', y='pct_exp_cells', by_exprs_values='counts', colour_by=ERCC_flag)+scale_x_log10()+
      ylab('percentage of expressing cells')+xlab('mean expression')+
      geom_hline(yintercept = 50, linetype='dashed')+geom_vline(xintercept = mean(assay(sce,'counts')),linetype='dashed')+
      geom_text(label=sprintf('%d genes are expressed in at least 50%% of cells', sum(rowData(sce)$pct_exp_cells>=50)), 
                y=40, x=10^(-05), size=6)+
      geom_text(label=sprintf('%d genes are expressed in at least 25%% of cells', sum(rowData(sce)$pct_exp_cells>=25)), 
                y=20, x=10^(-05), size=6)+
      geom_text(label=sprintf('%d genes are expressed in at least 1 cell', sum(rowData(sce)$detected>0)), 
                y=60, x=10^(-05), size=6)
    plot(plot_expr_vs_mean)
    dev.off()
  }
  # plot3: ggsave. mean expression vs frequency expressed: DOESN'T WORK
  if (F){
    pngfile=paste0("QC_meanexp_vs_freq",plate_id,".png")
    ggsave(pngfile, device='png')
    unlink(pngfile)
    # plot_expr_vs_mean <-  plotExprsFreqVsMean(sce, exprs_values = "counts", show_smooth = F)
    # function plotExprsFreqVsMean is discontinuous, so make our own
    plot_expr_vs_mean <- plotRowData(sce, x='mean', y='pct_exp_cells', by_exprs_values='counts', colour_by=ERCC_flag)+scale_x_log10()+
      ylab('percentage of expressing cells')+xlab('mean expression')+
      geom_hline(yintercept = 50, linetype='dashed')+
      geom_text(label=sprintf('%d genes are expressed in at least 50%% of cells', sum(rowData(sce)$pct_exp_cells>=50)), 
                y=40, x=10^(-05), size=6)+
      geom_text(label=sprintf('%d genes are expressed in at least 25%% of cells', sum(rowData(sce)$pct_exp_cells>=25)), 
                y=20, x=10^(-05), size=6)+
      geom_text(label=sprintf('%d genes are expressed in at least 1 cell', sum(rowData(sce)$detected>0)), 
                y=60, x=10^(-05), size=6)
    plot(plot_expr_vs_mean)
    png(pngfile, pointsize = 12, width = 720, height = 720)
    dev.off()
  }
  # ===== 1.2 QC_violinplots_X.pdf ====
  pdf(paste0("QC_violinplots_",plate_id,".pdf"), pointsize = 12)
  pviolin_idx <- 1
  plot_violin_scatter <- list()
  
  # plot for plates
  # a) sum
  plot_violin_scatter[[pviolin_idx]] <- plotColData(sce, x = "plate_id", y = "sum",
                                    colour_by = cellctrl_flag, size_by=nbcell_flag) +
    xlab("plate_id") +
    #stat_smooth(method = "lm", se = FALSE, size = 1.5, fullrange = TRUE) +
    ggtitle("total counts")
  if (cellsize_flip_flag) plot_violin_scatter[[pviolin_idx]] <- plot_violin_scatter[[pviolin_idx]] + scale_size_manual(values=nbcell_markersize)
  
  plot(plot_violin_scatter[[pviolin_idx]])
  pviolin_idx <- pviolin_idx+1
  
  # b) detected
  plot_violin_scatter[[pviolin_idx]] <- plotColData(sce, x = "plate_id", y = "detected",
                                                    colour_by = cellctrl_flag, size_by=nbcell_flag) +
    xlab("plate_id") +
    #stat_smooth(method = "lm", se = FALSE, size = 1.5, fullrange = TRUE) +
    ggtitle("detected transcripts")
  if (cellsize_flip_flag) plot_violin_scatter[[pviolin_idx]] <- plot_violin_scatter[[pviolin_idx]] + scale_size_manual(values=nbcell_markersize)
  plot(plot_violin_scatter[[pviolin_idx]])
  pviolin_idx <- pviolin_idx+1
  
  # c) pc_MT
  plot_violin_scatter[[pviolin_idx]] <- plotColData(sce, x = "plate_id", y = "subsets_mitochondrial_percent",
                                                    colour_by = cellctrl_flag, size_by=nbcell_flag) +
    xlab("plate_id") +
    #stat_smooth(method = "lm", se = FALSE, size = 1.5, fullrange = TRUE) +
    ggtitle("mitochondrial expression percentage ~ total feature")
  if (cellsize_flip_flag) plot_violin_scatter[[pviolin_idx]] <- plot_violin_scatter[[pviolin_idx]] + scale_size_manual(values=nbcell_markersize)
  plot(plot_violin_scatter[[pviolin_idx]])
  pviolin_idx <- pviolin_idx+1
  
  # c) pc_ERCC
  if (!is.null(ERCC_flag)){
    plot_violin_scatter[[pviolin_idx]] <- plotColData(sce, x = "plate_id", y = "subsets_ERCC_percent",
                                                      colour_by = cellctrl_flag, size_by=nbcell_flag) +
      xlab("plate_id") +
      #stat_smooth(method = "lm", se = FALSE, size = 1.5, fullrange = TRUE) +
      ggtitle("ERCC expression percentage ~ total feature")
    if (cellsize_flip_flag) plot_violin_scatter[[pviolin_idx]] <- plot_violin_scatter[[pviolin_idx]] + scale_size_manual(values=nbcell_markersize)
    plot(plot_violin_scatter[[pviolin_idx]])
    pviolin_idx <- pviolin_idx+1
  }
  
  
  # plot on meta1 / meta2 as category type
  for (meta_col in c('meta1','meta2')){
    if (meta_col=='meta1'){meta_flag=meta1_flag; xlab_txt='meta1'}
    if (meta_col=='meta2'){meta_flag=meta2_flag; xlab_txt='meta2'}
    if (!is.null(meta_flag)) {
      # a) sum
      plot_violin_scatter[[pviolin_idx]] <- plotColData(sce, x = xlab_txt, y = "sum",
                                                        colour_by = cellctrl_flag, size_by=nbcell_flag) +
        xlab(xlab_txt) +
        #stat_smooth(method = "lm", se = FALSE, size = 1.5, fullrange = TRUE) +
        ggtitle("total counts")
      if (cellsize_flip_flag) plot_violin_scatter[[pviolin_idx]] <- plot_violin_scatter[[pviolin_idx]] + scale_size_manual(values=nbcell_markersize)
      plot(plot_violin_scatter[[pviolin_idx]])
      pviolin_idx <- pviolin_idx+1
      
      # b) detected
      plot_violin_scatter[[pviolin_idx]] <- plotColData(sce, x = xlab_txt, y = "detected",
                                                        colour_by = cellctrl_flag, size_by=nbcell_flag) +
        xlab(xlab_txt) +
        #stat_smooth(method = "lm", se = FALSE, size = 1.5, fullrange = TRUE) +
        ggtitle("detected transcripts")
      if (cellsize_flip_flag) plot_violin_scatter[[pviolin_idx]] <- plot_violin_scatter[[pviolin_idx]] + scale_size_manual(values=nbcell_markersize)
      
      plot(plot_violin_scatter[[pviolin_idx]])
      pviolin_idx <- pviolin_idx+1
      
      # c) pc_MT
      plot_violin_scatter[[pviolin_idx]] <- plotColData(sce, x = xlab_txt, y = "subsets_mitochondrial_percent",
                                                        colour_by = cellctrl_flag, size_by=nbcell_flag) +
        xlab(xlab_txt) +
        #stat_smooth(method = "lm", se = FALSE, size = 1.5, fullrange = TRUE) +
        ggtitle("mitochondrial expression percentage ~ total feature")
      if (cellsize_flip_flag) plot_violin_scatter[[pviolin_idx]] <- plot_violin_scatter[[pviolin_idx]] + scale_size_manual(values=nbcell_markersize)
      
      plot(plot_violin_scatter[[pviolin_idx]])
      pviolin_idx <- pviolin_idx+1
      
      # c) pc_ERCC
      if (!is.null(ERCC_flag)){
        plot_violin_scatter[[pviolin_idx]] <- plotColData(sce, x = xlab_txt, y = "subsets_ERCC_percent",
                                                          colour_by = cellctrl_flag, size_by=nbcell_flag) +
          xlab(xlab_txt) +
          #stat_smooth(method = "lm", se = FALSE, size = 1.5, fullrange = TRUE) +
          ggtitle("ERCC expression percentage ~ total feature")
        if (cellsize_flip_flag) plot_violin_scatter[[pviolin_idx]] <- plot_violin_scatter[[pviolin_idx]] + scale_size_manual(values=nbcell_markersize)
        
        plot(plot_violin_scatter[[pviolin_idx]])
        pviolin_idx <- pviolin_idx+1
      }
    }
  }
  dev.off()
  # ===== 2. QC_dimReductions_X.pdf =====
  plot_pca <- list()
  pcaidx <- 1
  
  control_shape <- NULL
  if(length(inds_cellctrl)) control_shape <- "control"
  control_scale <- NULL
  if(length(inds_cellctrl)) control_scale <- scale_shape_manual(values=c(16,9))
  
  sce <- logNormCounts(sce)
  # pca on lognorm exp by default anyway
  sce <- BiocSingular::runPCA(sce, ncomponents=20, exprs_values = "logcounts")
  # UMAP
  sce <- runUMAP(sce, dimred = "PCA", ncomponents = 2, exprs_values = "logcounts",
                 name = "UMAP_on_PCA")
  # old version using t-SNE 
  if (F){
    set.seed(1000)
    sce <- scater::runTSNE(sce, perplexity=10)
    if (is_test) str(reducedDim(sce,'PCA'))
    
    pdf(paste0("QC_dimReductions_",plate_id,".pdf"), pointsize=12)
    
    plot_pca[[pcaidx]] <- plotPCA(sce,colour_by = "sum",size_by = "detected",
                        shape_by=control_shape) +control_scale+
      ggtitle('PCA on expression matrix')
    # plot_pca <- plot_pca + scale_fill_continuous(name='total counts')+
    #   scale_size_continuous(name = "total features")
    plot(plot_pca[[pcaidx]])
    pcaidx <- pcaidx+1
    # TSNE: can run multiple times with different seed and perplexity.
    
    plot_tsne <- plotTSNE(sce,
                          colour_by = "sum",size_by = "detected",shape_by = control_shape) + 
      control_scale+ggtitle('TSNE on expression matrix')
    plot(plot_tsne)
    
    # plot MT%, control%, ERCC% by TSNE if any
    plot_pca_ctr <- list()
    pidx <- 1
    plot_pca_ctr[[pidx]] <- plotPCA(sce,
                          colour_by = "subsets_mitochondrial_percent",size_by = "detected",shape_by = control_shape) + 
      control_scale+ggtitle('PCA on expression matrix, %MT per cell')
    plot(plot_pca_ctr[[pidx]])
    pidx <- pidx+1
    
    plot_pca_ctr[[pidx]] <- plotPCA(sce,
                          colour_by = "subsets_ERCC_percent",size_by = "detected",shape_by = control_shape) + 
      control_scale+ggtitle('PCA on expression matrix, %ERCC per cell')
    plot(plot_pca_ctr[[pidx]])
    pidx <- pidx+1
    dev.off()
  }
  
  pdf(paste0("QC_dimReductions_",plate_id,".pdf"), pointsize=12)
  # PCA
  plot_pca[[pcaidx]] <- plotPCA(sce,colour_by = "sum",size_by = "detected",
                                shape_by=control_shape) +control_scale+
    ggtitle('PCA on expression matrix')
  # plot_pca <- plot_pca + scale_fill_continuous(name='total counts')+
  #   scale_size_continuous(name = "total features")
  plot(plot_pca[[pcaidx]])
  pcaidx <- pcaidx+1
  
  # UMAP
  plot_pca[[pcaidx]] <- plotReducedDim(sce,'UMAP_on_PCA', colour_by = "sum",size_by = "detected",
                                shape_by=control_shape, ncomponents = 2) +control_scale+
    ggtitle('UMAP on expression matrix')
  plot(plot_pca[[pcaidx]])
  pcaidx <- pcaidx+1
  
  # plot MT%, control%, ERCC% by UMAP if any
  plot_pca_ctr <- list()
  pidx <- 1
  
  # if number of wells > 500, it becomes difficult to visualise with different marker size
  tmp.size_by <- NULL
  if (length(sce$unique_sample_id_suffix)<500) tmp.size_by <- 'detected'
  
  plot_pca_ctr[[pidx]] <- plotReducedDim(sce,'UMAP_on_PCA',colour_by = "subsets_mitochondrial_percent",
                                         size_by = tmp.size_by,shape_by = control_shape) + 
    control_scale+ggtitle('UMAP on expression matrix, %MT per cell')
  plot(plot_pca_ctr[[pidx]])
  pidx <- pidx+1
  
  plot_pca_ctr[[pidx]] <- plotReducedDim(sce,'UMAP_on_PCA',colour_by = "subsets_ERCC_percent",
                    size_by = tmp.size_by,shape_by = control_shape) + 
    control_scale+ggtitle('UMAP on expression matrix, %ERCC per cell')
  plot(plot_pca_ctr[[pidx]])
  pidx <- pidx+1
  
  dev.off()
  
  # ===== 3. QC_dimReductions_X.meta.pdf =====
  # old version with PCA
  if (F){
    if (samplesheet.check$meta_1 | samplesheet.check$meta_2){
      pdf(paste0("QC_dimReductions_",plate_id,".meta.pdf"), pointsize=12)
      colour_by=NULL
      if(samplesheet.check$meta_1 & !samplesheet.check$meta_2) {
        plot_pca[[pcaidx]] <- plotPCA(sce, colour_by = "sum", size_by = "detected", shape_by='meta1') +
          ggtitle('PCA on expression matrix using meta group for shape')
        plot(plot_pca[[pcaidx]])
        pcaidx <- pcaidx+1
      } else if (!samplesheet.check$meta_1 & samplesheet.check$meta_2) {
        plot_pca[[pcaidx]] <- plotPCA(sce, colour_by = "sum", size_by = "detected", shape_by='meta2') +
          ggtitle('PCA on expression matrix using meta group for shape')
        plot(plot_pca[[pcaidx]])
        pcaidx <- pcaidx+1
      } else if (samplesheet.check$meta_1 & samplesheet.check$meta_2) {
        plot_pca[[pcaidx]] <- plotPCA(sce, colour_by = "meta1", size_by = "detected", shape_by='meta2') +
          ggtitle('PCA on expression matrix using meta groups for colour and shape')
        plot(plot_pca[[pcaidx]])
        pcaidx <- pcaidx+1
      }
      dev.off()
    }
  }
  if (samplesheet.check$meta_1 | samplesheet.check$meta_2){
    tmp.size_by <- NULL
    if (length(sce$unique_sample_id_suffix)<500) tmp.size_by <- 'detected'
    
    pdf(paste0("QC_dimReductions_",plate_id,".meta.pdf"), pointsize=12)
    colour_by=NULL
    if(samplesheet.check$meta_1 & !samplesheet.check$meta_2) {
      plot_pca[[pcaidx]] <- plotReducedDim(sce,'UMAP_on_PCA', colour_by = "sum", size_by = tmp.size_by, shape_by='meta1') +
        ggtitle('UMAP on expression matrix using meta group for shape')
      plot(plot_pca[[pcaidx]])
      pcaidx <- pcaidx+1
    } else if (!samplesheet.check$meta_1 & samplesheet.check$meta_2) {
      plot_pca[[pcaidx]] <- plotReducedDim(sce,'UMAP_on_PCA', colour_by = "sum", size_by = tmp.size_by, shape_by='meta2') +
        ggtitle('UMAP on expression matrix using meta group for shape')
      plot(plot_pca[[pcaidx]])
      pcaidx <- pcaidx+1
    } else if (samplesheet.check$meta_1 & samplesheet.check$meta_2) {
      plot_pca[[pcaidx]] <- plotReducedDim(sce,'UMAP_on_PCA', colour_by = "meta1", size_by = tmp.size_by, shape_by='meta2') +
        ggtitle('UMAP on expression matrix using meta groups for colour and shape')
      plot(plot_pca[[pcaidx]])
      pcaidx <- pcaidx+1
    }
    dev.off()
    
  }
  
  # ===== 4. QC_dimReductions_X.plate_id.pdf =====
  if(length(unique(sce$plate_id))>1) {
    pdf(paste0("QC_dimReductions_",plate_id,"plate_id.pdf"), pointsize=12)
    # control_scale <- scale_shape_manual(values=c(16,9))
    if (length(unique(sce$plate_id))<10){
      plot_pca[[pcaidx]] <- plotReducedDim(sce,'UMAP_on_PCA',
                                           colour_by = "sum",
                                           size_by = "detected",
                                           shape_by='plate_id')  #+ control_scale
        
    } else{
      plot_pca[[pcaidx]] <- plotReducedDim(sce,'UMAP_on_PCA',
                                           colour_by = "plate_id")  #+ control_scale
    }
    plot_pca[[pcaidx]] <- plot_pca[[pcaidx]]+ggtitle('UMAP on expression matrix using plate id as shape')
    plot(plot_pca[[pcaidx]])
    pcaidx <- pcaidx+1
    dev.off()
  }
  
  # ===== 5. Plate_position_plots_X.pdf ====
  pdf(paste0("Plate_position_plots_",plate_id,".pdf"), pointsize = 12)
  plot_plate_counts <- plotPlatePosition(sce, colour_by = "sum", by_exprs_values = "counts", point_size = 9) + 
    guides(fill="colorbar") + theme(legend.position = "top")#+ggtitle('total counts on each well')
  plot(plot_plate_counts)
  plot_plate_mito <- plotPlatePosition(sce, colour_by = "subsets_mitochondrial_percent", by_exprs_values = "counts", point_size = 9) + 
    guides(fill="colorbar") + theme(legend.position = "top")#+ggtitle('%MT count on each well')
  plot(plot_plate_mito)
  plot_plate_ercc <- plotPlatePosition(sce, colour_by = "subsets_ERCC_percent", by_exprs_values = "counts", point_size = 9) + 
    guides(fill="colorbar")+ theme(legend.position = "top")#+ggtitle('%ERCC count on each well')
  plot(plot_plate_ercc)
  dev.off()
  # ===== 6. Sample_counts_across_ERCCs_X.tsv =====
  # ERCC mol
  # >add %ERCC
  ercc_counts <- counts(sce)[inds_ercc,] %>% as.data.frame(stringsAsFactors=F)
  
  # colnames(ercc_counts) <- ercc_counts %>% colnames() %>% gsub("_.*$","",.)
  ercc_counts$All_samples <- rowSums(ercc_counts)
  ercc_counts$ERCC <- rownames(ercc_counts)
  df_col <- ncol(ercc_counts)
  
  ercc_table <- ercc_counts[,c(df_col, df_col-1, 1:(df_col-2))]
  colnames(ercc_table) <- c('ERCC', 'All_samples',sce$fastqname)
  
  write.table(ercc_table, paste0("Sample_counts_across_ERCCs_",plate_id,".tsv"),
              sep='\t', row.names = F)

  # ==== 7. outlier detection ====
  
  # outlier cell condition 1: MT >10%
  MT_thresh <- 10
  sce$outlier_pct_MT <- sce$subsets_mitochondrial_percent>MT_thresh
  
  # condition 2: read mapping rate < 50%
  # read alignment rate file
  maprate_df <- read.table(pct_pseudo_txt, header = F, check.names = F)
  colnames(maprate_df) <- c('readname','map_rate')
  maprate_df <- maprate_df[-(maprate_df$readname=='mean'),]  
  # change readname eg: CU5DAY0_R0759-S0001_A64199_CU5DAY0A10_HTKKKDRXY_GGAGCTAC-AGAGGATA_L002 -> CU5DAY0A10
  if (F) {
    replace_readname <- function(newname, oldnames){
      return(which(grepl(pattern = paste0(newname,'_'), oldnames) ))
      # return index to oldname element: oldnames[idx] <- newname
    }
    idx_oldnames <- sapply(X = colnames(sce), FUN = replace_readname, oldnames=maprate_df$readname)
    maprate_df$Cell <- ""
    maprate_df$Cell[idx_oldnames] <- colnames(sce)
  }
  
  # add mapping rate to sce
  sce$map_rate <- 0
  sce$map_rate <- sapply(sce$fastqname, 
                         FUN = function(colname) maprate_df$map_rate[which(colname==maprate_df$readname)])
  sce$outlier_mapping_rate <- sce$map_rate<50
  
  sce$outlier_hits <- sce$outlier_pct_MT + sce$outlier_mapping_rate
  
  #saving per cell stats in a table
  per_cell_metric <- c("plate_id","plate_position", "outlier_hits","fastqname",
                       "outlier_pct_MT" ,"outlier_mapping_rate", "sum", "detected","unique_sample_id_suffix",
                       "subsets_mitochondrial_percent", "subsets_ERCC_percent", "map_rate")
  
  per_cell_table <- as.data.frame(colData(sce)[, per_cell_metric])
  per_cell_table$is_cell_control = F
  per_cell_table$is_cell_control[inds_cellctrl] <- T
  
  if (F) { # Vlad's setting
    rownames(per_cell_table) <- rownames(per_cell_table) %>% 
      gsub("_[CATG]{8}-[CATG]{8}","",.)
    sample_names <- per_cell_table %>% rownames
    per_cell_table$sample_ID <- sample_names
    per_cell_table <- per_cell_table %>% select(sample_ID, everything())
    write.table(per_cell_table, paste0("Per_sample_key_metrics_",plate_id,".tsv"),
                sep='\t', row.names = F)
    
    double_outliers <- which(per_cell_table$outlier_hits==2) %>% sample_names[.]
    single_outliers <- which(per_cell_table$outlier_hits==1) %>% sample_names[.]
  }
  
  rownames(per_cell_table) <- per_cell_table$fastqname
  double_outliers <- which(per_cell_table$outlier_hits==2) %>% per_cell_table$unique_sample_id_suffix[.]
  single_outliers <- which(per_cell_table$outlier_hits==1) %>% per_cell_table$unique_sample_id_suffix[.]
  per_cell_table <- subset(per_cell_table, select=-c(unique_sample_id_suffix))
  write.table(per_cell_table, paste0("Per_sample_key_metrics_",plate_id,".tsv"),
              sep='\t', row.names = F)
  
  # plot Mapping_rate_X.pdf
  pdf(paste0("Mapping_rate_",plate_id,".pdf"), pointsize = 12)
  cellctrl_flag <- NULL
  if(length(inds_cellctrl)>0){
    cellctrl_flag <- "control"
  }
  
  meta1_flag <- NULL
  if (length(unique(sce$meta1))>1) meta1_flag <- 'meta1'
  meta2_flag <- NULL
  if (length(unique(sce$meta2))>1) meta2_flag <- 'meta2'
  
  # optimal display:
  call_plotColData <- function(sce=NULL, x = "sum",y = "map_rate",
                               colour_by = NULL, size_by= NULL, shape_by=NULL, 
                               cellsize_flip_flag=F, nbcell_markersize=NULL) {
    
    p <- plotColData(sce, x = x,y = y,
                                                  colour_by = colour_by, size_by= size_by, shape_by=shape_by)
    p <- p+scale_x_log10()
    if (cellsize_flip_flag) p <- p + scale_size_manual(values=nbcell_markersize)
    
    p <- p+theme(legend.position = "right") + xlab("log_total_counts") + ylab("kallisto mapping rate") +
      labs(title="Read mapping rate ~ total counts (log scale)")
    return(p)
  }
  plot_mapping_rate <- list()
  idx_plot <- 1
  if (!is.null(cellctrl_flag)){
    colour_by = cellctrl_flag
    size_by= nbcell_flag
    shape_by=NULL
    if (!is.null(meta1_flag)){
      shape_by=meta1_flag
      plot_mapping_rate[[idx_plot]] <- call_plotColData(sce, colour_by = colour_by, size_by= size_by,
                                                        cellsize_flip_flag=cellsize_flip_flag, shape_by=meta1_flag, nbcell_markersize=nbcell_markersize)
      plot(plot_mapping_rate[[idx_plot]])
      idx_plot <- idx_plot+1
    }
    if (!is.null(meta2_flag)){
      shape_by=meta2_flag
      plot_mapping_rate[[idx_plot]] <- call_plotColData(sce, colour_by = colour_by, size_by= size_by,
                                                        cellsize_flip_flag=cellsize_flip_flag, shape_by=meta2_flag,nbcell_markersize=nbcell_markersize)
      plot(plot_mapping_rate[[idx_plot]])
      idx_plot <- idx_plot+1
    }
    if(is.null(meta1_flag) &  is.null(meta2_flag)){
      plot_mapping_rate[[idx_plot]] <- call_plotColData(sce, colour_by = colour_by, size_by= size_by,
                                                        cellsize_flip_flag=cellsize_flip_flag, shape_by=NULL, nbcell_markersize=nbcell_markersize)
      plot(plot_mapping_rate[[idx_plot]])
      idx_plot <- idx_plot+1
    }
  } else {
    size_by=nbcell_flag
    shape_by=NULL
    colour_by=NULL
    if (!is.null(meta1_flag)){
      plot_mapping_rate[[idx_plot]] <- call_plotColData(sce, colour_by=meta1_flag, size_by= size_by,
                                                        cellsize_flip_flag=cellsize_flip_flag, shape_by=shape_by, nbcell_markersize=nbcell_markersize)
      plot(plot_mapping_rate[[idx_plot]])
      idx_plot <- idx_plot+1
      }
    if (!is.null(meta2_flag)){
      plot_mapping_rate[[idx_plot]] <- call_plotColData(sce, colour_by=meta2_flag, size_by= size_by,
                                                        cellsize_flip_flag=cellsize_flip_flag, shape_by=shape_by, nbcell_markersize=nbcell_markersize)
      plot(plot_mapping_rate[[idx_plot]])
      idx_plot <- idx_plot+1
    }
    if(is.null(meta1_flag) &  is.null(meta2_flag)){
      plot_mapping_rate[[idx_plot]] <- call_plotColData(sce, colour_by=NULL, size_by= size_by,
                                                        cellsize_flip_flag=cellsize_flip_flag, shape_by=shape_by, nbcell_markersize=nbcell_markersize)
      plot(plot_mapping_rate[[idx_plot]])
      idx_plot <- idx_plot+1
    }
  }
  dev.off()

  
# ===== 8. save work space for doc compilation ====
  save(plot_cumulative_dist, plot_highest_expression, plot_expr_vs_mean, plot_featurenum_vs_counts, 
       plot_ercc_scatter, plot_mito_scatter,
       plot_plate_counts, plot_plate_mito, plot_plate_ercc,plot_violin_scatter,plot_mapping_rate,
       plot_pca, plot_pca_ctr,
       per_cell_table, ercc_table, per_cell_metric,
       double_outliers, single_outliers,
       plate_id,inds_cellctrl,
       sce,
       file=paste0(plate_id,"qc_for_doc.Rdata"))
  
  metadata.cell <- colData(sce)
  rownames(metadata.cell) <- metadata.cell$fastqname;
  write.table(x = metadata.cell, sep = '\t',
              file = paste0(plate_id, ".metadata_cell.tsv"))
  
  metadata.transcript <- rowData(sce)
  write.table(x = metadata.transcript, sep = '\t',
              file = paste0(plate_id, ".metadata_transcript.tsv"))
  
  # test loading into Seurat
  if (F){
    # for human control
    output_location="/Volumes/core-bioinformatics/development/lany/CB-GENANNO-523_Karim_Gharbi_EI_KG_ENQ-3708/Analysis/scqc_reqs-0.2.2/run1/qc_dir/"
    plate_id <- "R0121-P0002"
    # load expression matrix
    raw_counts <- read.table(file=file.path(output_location, '../quants_dir',
                                            paste0("est_counts",plate_id,"_matrix.tsv")), 
                             header = T, check.names = F)
    # for cowseq3
    plate_id <- "all"
    output_location <- "/Volumes/core-bioinformatics/analysis/CB-GENANNO-525_Charlotte_Utting_EI_CU_ENQ-5286_A_01/Analysis/scqc_reqs-0.2.1/run1/qc_dir"
    # load expression matrix
    raw_counts <- read.table(file=file.path(output_location, '../quants_dir',
                                            "all_plates.tsv"), header = T, check.names = F)
    # for PPBFX-1067_Claudia_Ribeiro_de_Almeida
    plate_id <- "all"
    output_location <- "/ei/cb/analysis/CB-PPBFX-1067_Claudia_Ribeiro_de_Almeida_BI_CR_ENQ-5220_A_01_Smart-Seq2_QC/Analysis/scqc_reqs-0.2.2/run2/qc_dir"
    # load expression matrix
    raw_counts <- read.table(file=file.path(output_location, '../quants_dir',
                                            "all_plates.tsv"), header = T, check.names = F)
    
    # for GENANNO-543_Stuart_Rushworth_UOE_SR
    plate_id <- "all"
    output_location <- "/ei/cb/analysis/CB-GENANNO-543_Stuart_Rushworth_UOE_SR_ENQ-5327_A_01/Analysis/scqc_reqs-0.2.2/run1/qc_dir"
    # load expression matrix
    raw_counts <- read.table(file=file.path(output_location, '../quants_dir',
                                            "all_plates.tsv"), header = T, check.names = F)
    library("Seurat")
    
    raw_counts[1:3,1:2]
    mydata <- CreateSeuratObject(counts=raw_counts, min.cells = 3, min.features = 200, project = 'human')
    
    # load meta data
    celltsvfile <- file.path(output_location, paste0(plate_id, ".metadata_cell.tsv"))
    cellmeta <- read.table(file=celltsvfile, header = T, check.names = F)
    dim(cellmeta)
    cellmeta[1:2,1:2]
    colnames(cellmeta)
    nrow(cellmeta)
    ncol(mydata)
    mydata <- AddMetaData(mydata, metadata = cellmeta)
    # show metadata
    head(mydata[[]])
    head(mydata@meta.data)
    # or add a selected meta data field "dummy"
    mydata$dummy <- "dummy"
    head(mydata@meta.data$dummy)
    
    # ===== 1.3: add feature meta data ====
    featuretsvfile <- file.path(output_location, paste0(plate_id, ".metadata_transcript.tsv"))
    featmeta <- read.table(file=featuretsvfile, header = T, check.names = F)
    featmeta[1:2,1:2]
    # > featmeta[1:2,1:2]
    #                       detected subsets_cellctrl_mean
    # TpMuguga-02g00001-t1 0.2604167                     0
    # TpMuguga-02g00002-t1 0.0000000                     0
    
    # keys for successfully loading cell meta data:
    # 1. featmeta is a dataframe, using feature as rownames, 
    #     some features may have been filtered out of the count matrix (mydata[[]])
    # 2. featmeta colnames are meta data fields.
    
    colnames(featmeta)
    # change TpMuguga_02g00003_t1 to TpMuguga-02g00003-t1
    newrowname <- gsub('_','-',rownames(featmeta)) %>% gsub('\\|','-',.)
    
    rownames(featmeta) <- newrowname
    featmeta[1:2,1:2]
    nrow(featmeta)
    nrow(mydata) # dim(mydata[["RNA"]])
    mydata[["RNA"]] <- AddMetaData(mydata[["RNA"]], metadata = featmeta)
    head(mydata[['RNA']][1:2,1:6])
    mydata[['RNA']]@meta.features[1:3,]
    
  }
  
  # test loading into seurat for PPBFX-1067 Claudia Ribeiro de Almeida BI.CR.ENQ-5220.A.01 Smart-Seq 2 QC
  if (F) {
    plate_id <- "all"
    output_location <- "/ei/cb/analysis/CB-PPBFX-1067_Claudia_Ribeiro_de_Almeida_BI_CR_ENQ-5220_A_01_Smart-Seq2_QC/Analysis/scqc_reqs-0.2.2/run2/qc_dir"
    # load expression matrix
    raw_counts <- read.table(file=file.path(output_location, '../quants_dir',
                                            "all_plates.tsv"), header = T, check.names = F)
    
    library("Seurat")
    
    raw_counts[1:3,1:2]
    mydata <- CreateSeuratObject(counts=raw_counts, min.cells = 3, min.features = 200, project = 'mouse')
    
    # load meta data
    celltsvfile <- file.path(output_location, paste0(plate_id, ".metadata_cell.tsv"))
    cellmeta <- read.table(file=celltsvfile, header = T, check.names = F)
    dim(cellmeta)
    cellmeta[1:2,1:2]
    colnames(cellmeta)
    nrow(cellmeta)
    ncol(mydata)
    mydata <- AddMetaData(mydata, metadata = cellmeta)
    # show metadata
    head(mydata[[]])
    head(mydata@meta.data)
    # or add a selected meta data field "dummy"
    mydata$dummy <- "dummy"
    head(mydata@meta.data$dummy)
    
    # ===== 1.3: add feature meta data ====
    featuretsvfile <- file.path(output_location, paste0(plate_id, ".metadata_transcript.tsv"))
    featmeta <- read.table(file=featuretsvfile, header = T, check.names = F)
    featmeta[1:2,1:2]
    # > featmeta[1:2,1:2]
    #                       detected subsets_cellctrl_mean
    # TpMuguga-02g00001-t1 0.2604167                     0
    # TpMuguga-02g00002-t1 0.0000000                     0
    
    # keys for successfully loading cell meta data:
    # 1. featmeta is a dataframe, using feature as rownames, 
    #     some features may have been filtered out of the count matrix (mydata[[]])
    # 2. featmeta colnames are meta data fields.
    
    colnames(featmeta)
    # change TpMuguga_02g00003_t1 to TpMuguga-02g00003-t1
    newrowname <- gsub('_','-',rownames(featmeta)) %>% gsub('\\|','-',.)
    
    rownames(featmeta) <- newrowname
    featmeta[1:2,1:2]
    nrow(featmeta)
    nrow(mydata) # dim(mydata[["RNA"]])
    mydata[["RNA"]] <- AddMetaData(mydata[["RNA"]], metadata = featmeta)
    head(mydata[['RNA']][1:2,1:6])
    mydata[['RNA']]@meta.features[1:3,]
    
  
}
}
if(sys.nframe()==0){
  outloc <- "./"
  sheetfile <- NULL
  mtrds <- NULL
  args=commandArgs(trailingOnly=T)
  
  if(!file.exists(args[1]))
    stop(paste0("Matrix file ('",args[1],"') not found!"))
  if(!file.exists(args[2]))
    stop(paste0("kallisto psudo align mapping file ('",args[2],"') not found!"))
  if(length(args)>=3){
    outloc <- args[3]
    if(!file.exists(outloc)){
      warning(paste0("Output directory ('",outloc,"') does not exist!"),
              "\nCreating the directory...")
      system(paste("mkdir -p",outloc))
    }
  }
  if(length(args)>=4)
    sheetfile <- args[4]
  if(length(args)>=5)
    mtrds <- args[5]

  scqc_from_tsv(args[1], args[2], outloc, sheetfile, mtrds)
}
