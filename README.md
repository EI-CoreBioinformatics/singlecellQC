# Smart-Seq2 quality control pipeline

Smart-seq2 quantification and quality control pipeline, developed for single-cell service at the Earlham Institue.
version 1.1

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

* kallisto 0.48.0
* nextflow
* R and R packages (scater 1.28.0, rjson 0.2.21, rmarkdown, tinytex, knitr, dplyr)
* singularity 2.4.2
* pandoc, texlive, ghostscript

## Preprocessing
1. convert the sample sheet to UNIX format
```
source dos2unix-7.4.1_CBG
dos2unix -n old_samplesheet samplesheet
```
2. format read names. The pipeline can take symbolic links and requires read names in the following format: 
```
${Sample_Plate}_${demultiplexed_readname}
``` 
For example, the following FASTQ file from plate `CUB35DAY0` with a filename:
```
R0882-S0001_A68701_CUB35DAY0A10_H3VY5DRX2_CGTACTAG-AGAGGATA_L001_R1.fastq.gz
```
becomes
```
CUB35DAY0_R0882-S0001_A68701_CUB35DAY0A10_H3VY5DRX2_CGTACTAG-AGAGGATA_L001_R1.fastq.gz
```

## Running the pipeline

Pipeline is written in Nextflow, so a run is usually initiated in the following way:
`nextflow run scqc_nf.sh -c config_file &`

Examples of a config file and sample sheet are in the repository.

An example:

```
cd /ei/cb/development/lany/CB-GENANNO-525_Charlotte_Utting_EI_CU_ENQ-5286_A_01/Analysis/scqc_reqs-1.1/4plates.run2
sbatch -p ei-cb -J scqc_GENANNO-525.all_plate -o scqc_GENANNO-525.all_plate.%j.%N.log -c 1 --mem 10G \
    --mail-type=ALL --mail-user=user.email.com \
    --wrap "source singlecellQC-1.1_CBG && cd $analysis_dir && \
    nextflow run /ei/software/cb/singlecellQC/1.1/x86_64/bin/scqc_nf.sh \
    -c GENANNO-525.scqc-1.1.all_plates.config -with-report -resume"
```

## Config file

Parameters inside 'params' scope can be passed when starting the pipeline by adding them at the end of the pipeline start call, e.g.
'nextflow run scqc_nf.sh -c scqc.config --qcoutdir=my_qc_directory'. Alternatively, they can be edited in the config file.

Parameters related to the output, organism species and this pipeline.

* quantificationsdir - Directory to contain qunatifications produced for each of the samples and the counts matrices
* qcoutdir' - Directory to contain the final QC report and other QC-related files
* samplesheet - .csv file containing information about sample names, wells, control status and other sample metadata.
* reads -  Location of sample FASTQ files.
* species - E.g., 'Hsapiens', 'Mmusculus'.
* plate_ids - List of plate identificators (typically 4 strings), as they appear in the names of raw data samples. This is how the pipeline merges
    samples into plate-level matrices which are then used for plate-level QC.
* mtnamefile - In case of non-human species, .rds file for mitochondrial gene. Leave empty ('') if human.
    This is a vector in R, containing the list of Ensembl transcript IDs, saved as .rds (using saveRDS())
* pattern - The format of the FASTQ endings showing how they should be grouped, as a glob pattern
* trans2gen_tsv - Transcript id to gene id mapping tsv file

General HPC parameters in scope `executor`:

* executor - Type of HPC scheduler.

* queue - Queue used for submitting jobs.

* memory - Memory assigned to a job.

* queueSize - Maximum number of jobs the pipeline will submit at once.

Process specific parameters can be set in scope `process`.

For more information on configuration file options, check out [nextflow documentation](https://www.nextflow.io/docs/latest/config.html).


## Sample sheet

Sample sheet will be unique for every run.

It is a .csv file that has to have the following columns. Additional columns are not a problem, but are not used. 
Make sure there is no white space in any entries.

* Sample_ID, Sample_Name - These two columns are also in the Illumina sample sheet.
* Sample_Plate - A string corresponding to one of the plates used in the experiment
* Sample_Well - Row/column location of the sample on the plate. Something like A01, A02, etc. These are needed for plate position plots.
* row - Row location of the sample on the plate. Normally ranges between A to H.
* column - Column location of the sample on the plate. Normal range is between 1 and 12.
* control - TRUE if the well contains a control, FALSE otherwise.
* number_of_cells - Number of cells in a well, e.g., 0, 1, 2, 20, 50. Wells contains more or less than 1 cells will be labeled as "control".
* meta_1: meta data field 1, catagorical value.
* meta_2: meta data field 2, catagorical value

## Required precomputed resources

* kallisto index
* list of mitochondrial genes
* transcript id to gene id mapping tsv file
* singularity images

