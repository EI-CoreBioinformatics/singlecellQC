# Smart-Seq2 quality control pipeline

Smart-seq2 quantification and quality control pipeline, developed for single-cell service at the Earlham Institue.

Starting with the FASTQ files generated from the SmartSeq2 experiment and a sample sheet containing the appropriate metadata, transcript level counts are generated, merged and used to produce QC metrics.
A QC report document is generated for each of the plates, and the whole experiment.

The pipeline does not require internet connection, but some required files have to be precomputed.

## Inputs

* Config file

## Outputs

* Sample count quantifications
* Expression counts matrix (transcript and gene-level)
* Quality control metrics and plots (per plate)
* QC report (per plate and for the entire experiment)


## Tools used

* kallisto 0.44.0
* nextflow 22.04.0
* R 4.1.2 and R packages (scater 1.22.0, ggplot2, dplyr, knitr, rjson)
* singularity 3.8.7


## Running the pipeline

Pipeline is written in Nextflow, so a run is usually initiated in the following way:
`nextflow run scqc_nf.sh -c config_file &`

As a batch job to HPC:
```
sbatch -p ei-cb -J jobname -o jobname.%j.log -c 1 --mem 10G \
    --mail-type=ALL --mail-user=user@email.com \
    --wrap "source /ei/cb/common/Scripts/singlecellQC/v1.0/scqc_reqs-1.0 && \
    nextflow run scqc_nf.sh -c config_file -with-report -resume"
```
Examples of a config file and sample sheet are in the repository.

## Config file

Parameters with 'params.' prefix can be passed when starting the pipeline by adding them at the end of the pipeline start call, e.g.
'nextflow run scqc_nf.sh -c scqc.config --qcoutdir=my_qc_directory'. Alternatively, they can be edited in the config file.

Parameters related to the output, organism species and this pipeline.

* quantificationsdir - Directory to contain qunatifications produced for each of the samples and the counts matrices

* qcoutdir - Directory to contain the final QC report and other QC-related files

* reads -  Location of sample FASTQ files.

* samplesheet - .csv file containing information about sample names, wells, control status and other sample metadata.

* plate_ids - List of plate identificators (typically 4 strings), as they appear in the names of raw data samples. This is how the pipeline merges
    samples into plate-level matrices which are then used for plate-level QC.
    
* species - Kept for legacy code, this parameter is only used if is 'Hsapiens' or 'Bos_taurus'.

* idx - Location of the kallisto index to use. Indices are precomputed. If you want to use a new one, one can be build with 'kallisto index'.

* biomartdb - Legacy variable name. This is a Ensembl transcript name to gene name tsv lookup table, and is required to convert expression count matrix         from transcript to gene level.

* mtnamefile - This is a vector in R, containing the list of Ensembl metochondrial transcript IDs, saved as .rds (using saveRDS()).  
    Leave empty ("") there is no MT annotation.
    
* pattern - The format of the FASTQ endings showing how they should be grouped, as a glob pattern

General HPC parameters:

* executor - Type of HPC scheduler.

* queue - Queue used for submitting jobs.

* memory - Memory assigned to a job.

* queueSize - Maximum number of jobs the pipeline will submit at once.

For more information on configuration file options, check out [nextflow documentation](https://www.nextflow.io/docs/latest/config.html).


## Sample sheet

Sample sheet will be unique for every run.

It is a .csv file that has to have the following columns. Additional columns are not a problem, but are not used. 
Make sure there is no white space in any entries.

* Sample_ID - As in Illumina sample sheet.
* Sample_Name - As in Illumina sample sheet.
* Sample_Plate - A customer supplied field, one of the plates used in the experiment. The string does not always correspond to Sample_ID.
* Sample_Well - Row/column location of the sample on the plate. Something like A01, A02, etc. These are needed for plate position plots.
* row - A customer supplied field, e.g., A, B, etc.
* column - A customer supplied field, e.g., 1, 2, etc.
* control - TRUE if the well contains a control, FALSE otherwise.
* number_of_cells - Number of cells in a well, e.g., 0, 1, 10, 20, 50.
* meta_1 - A customer supplied meta data, categorical value.  
* meta_2 - A customer supplied meta data, categorical value. 

## Required precomputed resources

* kallisto index
* list of mitochondrial genes
* transcript to gene lookup table
* singularity images

