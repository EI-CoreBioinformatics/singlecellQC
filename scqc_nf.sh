// before running:'source scqc_reqs-0.1'
// standard run case:'nextflow run some_path/scqc_nf.sh -c config_file &'

//human_index='/ei/cb/common/Scripts/scqc/References/human_w_ercc.idx'
//mouse_index='/ei/cb/common/Scripts/scqc/References/mm_ERCC.idx'
//chicken_index='/ei/cb/common/Scripts/scqc/References/chicken_ERCC.idx'

image_scater='/ei/cb/common/Scripts/scqc/Containers/R-3.5.2_scater.img'
image_json='/ei/cb/common/Scripts/scqc/Containers/R-3.5.2_bioMjsonS3.img'
image_doc='/ei/cb/common/Scripts/scqc/Containers/verse.simg'

qc_script='/ei/cb/common/Scripts/scqc/Scripts/scqc_from_matrix.R'
merge_script='/ei/cb/common/Scripts/scqc/Scripts/merge_kallisto_quant.R'
k_scrape_script='/ei/cb/common/Scripts/scqc/Scripts/kallisto_mapping_scrape.R'
doc_script='/ei/cb/common/Scripts/scqc/Scripts/QCreport.Rmd'
plate_matrix_merge_script='/ei/cb/common/Scripts/scqc/Scripts/plate_merge.R'
tx2g_script='/ei/cb/common/Scripts/scqc/Scripts/est_counts_tx2gene.R'

k_scrape_metric='p_pseudoaligned'
k_scrape_output='percent_pseudoaligned.txt'
doc_rdata='qc_for_doc.Rdata'

Channel
	.fromFilePairs("${params.reads}/${params.pattern}.fastq.gz")
	.set{read_pairs}

process quantification {

    beforeScript 'source kallisto-0.45.1'
    publishDir "$params.quantificationsoutdir"
    tag "$sampleId"

    input:
    set sampleId, file(reads) from read_pairs
	
    output:
    file "${sampleId}" into quants

    """
    kallisto quant -i ${params.idx} -o $sampleId -b 100 $reads
    """
}

process p_kal {

	publishDir "$params.qcoutdir", mode: 'copy'
	
	input:
	file flag_check from quants.collect()

	output:
	file "${k_scrape_output}" into k_scrape_ch

	"""
	singularity exec ${image_json} Rscript ${k_scrape_script} ${params.quantificationsoutdir} ${k_scrape_output} ${k_scrape_metric}
	"""
}

process q_merge {
	tag "$plate_id"	
	errorStrategy 'finish'
	beforeScript 'source R-3.5.2'	
	
	input:
	file flag_check from k_scrape_ch
	val plate_id from params.plate_ids
	
	output:
	file "tsv_name_${plate_id}" into count_file

	"""	
	Rscript $merge_script ${params.quantificationsoutdir} est_counts ${plate_id};
	echo \$(ls -d -1 ${params.quantificationsoutdir}*.* | grep ${plate_id}_matrix.tsv) > tsv_name_${plate_id};
	"""
}

process qc {
	errorStrategy 'finish'
	beforeScript 'export HDF5_DISABLE_VERSION_CHECK=1'

	input:
	file name_file from count_file

	output:
	file name_file into qc_done

	""" 
	singularity exec ${image_scater} Rscript ${qc_script} \$(cat ${name_file}) ${params.qcoutdir} ${params.plate_info} ${params.samplesheet} ${params.mtnamefile}
	"""
}

process doc {
	publishDir "$params.qcoutdir", mode: 'copy'

	input:
	file flag2_check from qc_done.collect()
	val plate_id from params.plate_ids
	
	
	output:
	file "Finished_${plate_id}.txt" into finished_ch

	"""
	singularity exec ${image_doc} Rscript -e \"options(warn=-1);objects<-\'${params.qcoutdir}${plate_id}${doc_rdata}\';mapping_file <- read.table(\'${params.qcoutdir}${k_scrape_output}\');rmarkdown::render(\'${doc_script}\', 'pdf_document', output_file=\'${plate_id}_QC_report.pdf\', output_dir=\'${params.qcoutdir}${plate_id}\')\";
	echo '${plate_id}' > Finished_${plate_id}.txt
	"""
}

process mat_merge {
	publishDir "$params.quantificationsoutdir", mode: 'copy'

	input:
	val id_list from finished_ch.collect()

	output:
	file 'matrix_location.txt' into fortx2g_ch
	
	"""
	singularity exec ${image_scater} Rscript ${plate_matrix_merge_script} ${params.quantificationsoutdir} \'${id_list}\';
	echo \$(ls -d -1 ${params.quantificationsoutdir}*.* | grep all) > matrix_location.txt
	"""

}

process tx2g {
	publishDir "$params.quantificationsoutdir", mode: 'copy'

	input:
	file tx_matrix_location from fortx2g_ch

	output:
	file tx_matrix_location into forAllQC_ch

	"""
	singularity exec ${image_json} Rscript ${tx2g_script} \$(cat ${tx_matrix_location}) ${params.quantificationsoutdir} ${params.species}
	"""

}

process all_qc {
	beforeScript 'export HDF5_DISABLE_VERSION_CHECK=1'

	input:
	file name_file from forAllQC_ch

	output:
	file 'qfolder_name' into all_folder_name

	""" 
	singularity exec ${image_scater} Rscript ${qc_script} \$(cat ${name_file}) ${params.qcoutdir} ${params.plate_info} ${params.samplesheet} ${params.mtnamefile};
	echo ${params.quantificationsoutdir} > qfolder_name;
	"""
}

process all_doc {
	publishDir "$params.qcoutdir", mode: 'copy'

	input:
	file qfold_name from all_folder_name	
	
	output:
	file "Finished_all.txt" into all_finished_ch

	"""
	singularity exec ${image_doc} Rscript -e \"options(warn=-1);objects<-\'${params.qcoutdir}all${doc_rdata}\';mapping_file <- read.table(\'${params.qcoutdir}${k_scrape_output}\');rmarkdown::render(\'${doc_script}\', 'pdf_document', output_file=\'all_QC_report.pdf\', output_dir=\'${params.qcoutdir}all\')\";
	echo 'All' > Finished_all.txt
	"""
}
