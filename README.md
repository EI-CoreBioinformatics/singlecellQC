# Smart-Seq2 quality control pipeline

Smart-seq2 quantification and quality control pipeline, developed for single-cell service at the Earlham Institue.

Starting with the FASTQ files generated from the SmartSeq2 experiment and a sample sheet containing the appropriate metadata, transcript level counts are generated, merged and used to produce QC metrics.
A QC report document is generated for each of the plates, and the whole experiment.

The pipeline does not require internet connection, but some required files have to be precomputed.

## Inputs

* Demulitplexed read data (.fastq.gz)
* Sample sheet (.csv)
* Config file

## Outputs

* Sample count quantifications
* Expression counts matrix (transcript and gene-level)
* Quality control metrics and plots (per plate)
* QC report (per plate and for the entire experiment)


## Tools used

* kallisto 0.45.1
* nextflow 19.04
* R 3.5.2 and R packages (scater, rmarkdown, tidyverse)
* singularity 2.4.2


## Running the pipeline

Pipeline is written in Nextflow, so a run is usually initiated in the following way:
`nextflow run scqc_nf.sh -c config_file &`

Examples of a config file and sample sheet are in the repository.

## Config file

Parameters with 'params.' prefix can be passed when starting the pipeline by adding them at the end of the pipeline start call, e.g.
'nextflow run scqc_nf.sh -c scqc.config --qcoutdir=my_qc_directory'. Alternatively, they can be edited in the config file.

Parameters related to the output, organism species and this pipeline specifics.

* quantificationsdir - Directory to contain qunatifications produced for each of the samples and the counts matrices
* qcoutdir' - Directory to contain the final QC report and other QC-related files
* reads -  Location of sample FASTQ files.
* species - 'Hsapiens' or 'Mmusculus'
* samplesheet - .csv file containing information about sample names, wells, control status and other sample metadata.
* plate_ids - List of plate identificators (typically 4 strings), as they appear in the names of raw data samples. This is how the pipeline merges
    samples into plate-level matrices which are then used for plate-level QC.
* mtnamefile - In case of non-human species, .rds file for mitochondrial gene. Leave empty ('') if human.
    This is a vector in R, containing the list of Ensembl transcript IDs, saved as .rds (using saveRDS())

* idx - Location of the kallisto index to use. Indices are precomputed. If you want to use a new one, one can be build with 'kallisto index'.
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

* unique_sample_id_suffix - Part of the FASTQ file name that uniquely matches to one sample. This is how the pipeline connects the raw data with sample sheet information.
* well - Row/column location of the sample on the plate. Something like A01, A02, etc. These are needed for plate position plots.
* plate_id - A string corresponding to one of the plates used in the experiment.
* number_of_cells - Number of cells in a well.
* control - TRUE if the well contains a control, FALSE otherwise.
* experiment - Name of the experiment. Can be anything as long as it's the same througout the column.

## Required precomputed resources

* kallisto index
* list of mitochondrial genes
* singularity images
* biomaRt annotation

