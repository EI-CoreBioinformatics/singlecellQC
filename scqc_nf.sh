// link config file with '-c config_file'

human_index='/ei/workarea/group-pb/CB-PPBFX-751_Assessment_of_SmartSeq2_pool_IPO5217_PIP-2297/Scripts/human_w_ercc.idx'
mouse_index='/ei/workarea/group-pb/EICB_SingleCell/Reference/Added_ERCC/mm_ERCC.idx'
chicken_index='/ei/workarea/group-pb/EICB_SingleCell/Reference/Added_ERCC/chicken_ERCC.idx'
image_scater='/ei/workarea/group-pb/CB-PPBFX-751_Assessment_of_SmartSeq2_pool_IPO5217_PIP-2297/Software/R-3.5.2_scater.img'
image_json='/ei/workarea/group-pb/CB-PPBFX-751_Assessment_of_SmartSeq2_pool_IPO5217_PIP-2297/Software/R-3.5.2_bioMjsonS3.img'
image_doc='/hpc-home/uzun/containing/verse.simg'
qc_script='/hpc-home/uzun/Scripts/dev_scqc_from_matrix.R' //dev
merge_script='/ei/workarea/group-pb/CB-PPBFX-751_Assessment_of_SmartSeq2_pool_IPO5217_PIP-2297/Scripts/merge_kallisto_quant.R'
k_scrape_script='/ei/workarea/group-pb/CB-PPBFX-751_Assessment_of_SmartSeq2_pool_IPO5217_PIP-2297/Scripts/kallisto_mapping_scrape.R'
k_scrape_metric='p_pseudoaligned'
k_scrape_output='percent_pseudoaligned.txt'
doc_rdata='qc_for_doc.Rdata'
doc_script='/hpc-home/uzun/scr/QCreport.Rmd'

Channel
	.fromFilePairs("${params.reads}/${params.pattern}.fastq.gz")
	.set{read_pairs}

process quantification {

    beforeScript 'source kallisto-0.45.1'
    publishDir "$params.quantificationsoutdir"
    tag "$sampleID"

    input:
    set sampleId, file(reads) from read_pairs
	
    output:
    file "${sampleId}" into quants

    """
    kallisto quant -i ${params.idx} -o $sampleId -b 100 $reads
    """
}

process q_merge {
	tag 'Quant_merge'	
	
	beforeScript 'source R-3.5.2'	
	
	input:
	file flag_check from quants.collect()
	
	output:
	file 'tsv_name' into count_file

	"""	
	Rscript $merge_script ${params.quantificationsoutdir} est_counts;
	echo \$(ls -d -1 ${params.quantificationsoutdir}*.* | grep tsv) > tsv_name;
	"""
}

process qc {

	beforeScript 'export HDF5_DISABLE_VERSION_CHECK=1'

	input:
	file 'tsv_name' from count_file

	output:
	file 'qfolder_name' into folder_name

	"""
	singularity exec ${image_scater} Rscript ${qc_script} \$(cat tsv_name) ${params.qcoutdir} ${params.plate_info} ${params.samplesheet} ${params.mtnamefile};
	echo ${params.quantificationsoutdir} > qfolder_name;
	"""
}

process p_kal {

	publishDir "$params.qcoutdir", mode: 'copy'
	
	input:
	file 'qfolder_name' from folder_name

	output:
	file "${k_scrape_output}" into k_scrape_ch

	"""
	singularity exec ${image_json} Rscript ${k_scrape_script} \$(cat qfolder_name) ${k_scrape_output} ${k_scrape_metric}
	"""
}

process doc {
	publishDir "$params.qcoutdir", mode: 'copy'

	input:
	file scrap from k_scrape_ch //percentpseudo
	
	output:
	file "Finished.txt" into final_ch

	"""
	singularity exec ${image_doc} Rscript -e \"objects<-\'${params.qcoutdir}${doc_rdata}\';mapping_file <- read.table(\'${scrap}\');rmarkdown::render(\'${doc_script}\', 'pdf_document', output_dir=\'${params.qcoutdir}\')\";
	echo 'Pipeline run done.' > Finished.txt
	"""
}
