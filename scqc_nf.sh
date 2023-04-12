// before running: 'source scqc_reqs-1.0'
// standard run case:'nextflow run scqc_nf.sh -c PIP-3144.config -with-dag flowchart.png -with-report -resume'

nextflow.enable.dsl=2

image_scater='--overlay /ei/cb/common/Scripts/singlecellQC/Containers/R-scater.v7.uwat_overlay.img:ro /ei/cb/common/Scripts/singlecellQC/Containers/R_scater_1.22.img'
image_json='/ei/cb/common/Scripts/singlecellQC/Containers/R-3.5.2_bioMjsonS3.img'
image_doc='/ei/cb/common/Scripts/singlecellQC/Containers/R_verse.dplyr_1.0.8.img'
image_gs='/ei/cb/common/Scripts/singlecellQC/Containers/ghost_script.img'
image_kallisto='/ei/cb/common/Scripts/singlecellQC/Containers/kallisto_0.44.0.img'

qc_script='/ei/cb/development/lany/git/singlecellQC/main_dev_branch/Scripts/scqc_from_matrix.meta.R' 
merge_script='/ei/cb/development/lany/git/singlecellQC/main_dev_branch/Scripts/merge_kallisto_quant.R'
k_scrape_script='/ei/cb/development/lany/git/singlecellQC/main_dev_branch/Scripts/kallisto_mapping_scrape.R'
doc_script='/ei/cb/development/lany/git/singlecellQC/main_dev_branch/Scripts/QCreport.Rmd'
plate_matrix_merge_script='/ei/cb/development/lany/git/singlecellQC/main_dev_branch/Scripts/plate_merge.R'
tx2g_script='/ei/cb/development/lany/git/singlecellQC/main_dev_branch/Scripts/est_counts_tx2gene.R'

k_scrape_metric='p_pseudoaligned'
k_scrape_output='percent_pseudoaligned.txt'
doc_rdata='qc_for_doc.Rdata'
q_merge_col='est_counts'

process quantification {
	/* to config for a process only, e.g, :
	queue 'ei-cb'
	memory '5G'
	time '1h'
	*/

  publishDir "$params.quantificationsoutdir"
  tag "$sampleId"

  input:
  	tuple val(sampleId), path(reads)

  output:
  	path "${sampleId}" , emit: quants

  """
  singularity exec ${image_kallisto} kallisto quant -i ${params.idx} -o $sampleId -b 100 $reads
  """
}

process p_kal {

	// output percent_pseudoaligned.txt to qc_dir/ by copying
	publishDir "$params.qcoutdir", mode: 'copy'
	
	input:
		// wait for all kallisto processes to finish 
		path quants_check 

	output:
		path "${k_scrape_output}"

	"""
	singularity exec ${image_json} Rscript ${k_scrape_script} ${params.quantificationsoutdir} ${k_scrape_output} ${k_scrape_metric}
	"""
}


process q_merge {
	queue 'ei-cb'
	memory '5G'


	tag "$plate_id"	
	errorStrategy 'finish'
	
	input:
		// wait for p_kal to output percent_pseudoaligned.txt - will be ln to current process dir, although it is not used.
		path p_kal
		val plate_id 
	
	output:
		path "est_counts${plate_id}_matrix.tsv" , emit: estcounts

	// R script generate, e.g., quants_dir/est_countsCU5DAY0_matrix.tsv
	"""	
	singularity exec ${image_doc} Rscript $merge_script ${params.quantificationsoutdir} ${q_merge_col} ${plate_id};
	ln -s -v -f ${params.quantificationsoutdir}est_counts${plate_id}_matrix.tsv est_counts${plate_id}_matrix.tsv
	"""
}

process qc {
	errorStrategy 'finish'
	beforeScript 'export HDF5_DISABLE_VERSION_CHECK=1'

	input:
		// e.g., quants_dir/est_countsCU5DAY0_matrix.tsv
		path est_counts_file 
		path pc_pseudoalign_file

	output:
		path "qc_${est_counts_file}.done"

	""" 
	singularity exec ${image_scater} Rscript ${qc_script} ${est_counts_file} ${pc_pseudoalign_file} ${params.qcoutdir} ${params.samplesheet} ${params.mtnamefile}
	echo '${est_counts_file} done. ' > qc_${est_counts_file}.done
	"""
}

process gs {
	errorStrategy 'finish'
	publishDir "$params.qcoutdir", mode: 'copy'
 

	input:
		path qc_complete_check
		val plate_id

	output:
		path "QC_meanexp_vs_freq${plate_id}.png"

	"""	
	singularity exec ${image_gs} gs \
    -dNOPAUSE -dQUIET -dBATCH -sDEVICE=png16m -sOutputFile=QC_meanexp_vs_freq${plate_id}.png -r256 \
    ${params.qcoutdir}/QC_meanexp_vs_freq${plate_id}.pdf
  """
}

process doc {
	queue 'ei-cb'
	memory '10G'

	publishDir "$params.qcoutdir", mode: 'copy'
	//cache false

	input:
	path qc_complete_check
	path gs_png_file 
	val plate_id 
	
	
	output:
	path "Finished_${plate_id}.txt" 

	"""
	singularity exec ${image_doc} Rscript -e \"options(warn=-1);objects<-\'${params.qcoutdir}${plate_id}${doc_rdata}\';mapping_file <- read.table(\'${params.qcoutdir}${k_scrape_output}\');rmarkdown::render(\'${doc_script}\', 'pdf_document', output_file=\'${plate_id}_QC_report.pdf\', output_dir=\'${params.qcoutdir}${plate_id}\')\";
	echo '${plate_id}' > Finished_${plate_id}.txt
	"""
}

process mat_merge {
	publishDir "$params.quantificationsoutdir", mode: 'copy'
  //cache false

	input:
	// e.g., Finished_LTHSCO2.txt. 
	path doc_complete_check 

	output:
	path 'all_plates.tsv'
	
	
	// singularity exec R_verse.v5.img Rscript plate_merge.R quants_dir/ \
	//	'Finished_CU5DAY0.txt Finished_CU7DAY0.txt Finished_CU5DAY7.txt Finished_CU7DAY7.txt';
	"""
	singularity exec ${image_doc} Rscript ${plate_matrix_merge_script} ${params.quantificationsoutdir} ${q_merge_col} \'${doc_complete_check}\';
	ln -s -v -f ${params.quantificationsoutdir}all_plates.tsv all_plates.tsv
	"""
}

process tx2g {
	publishDir "$params.quantificationsoutdir", mode: 'copy'

	input:
	path tx_matrix_location 

	output:
	path 'tx2g_finished.txt'

	"""
	singularity exec ${image_json} Rscript ${tx2g_script} ${tx_matrix_location} ${params.quantificationsoutdir} ${params.species} ${params.biomartdb} 
	echo 'tx2g finished' > tx2g_finished.txt
	"""

}

workflow QC_AND_DOC {
	take: 
		est_counts_file 
		pc_pseudoalign_file
		plate_ids
	main:
		qc( est_counts_file, pc_pseudoalign_file )
	  gs(qc.out.collect(), plate_ids)
	  doc( qc.out.collect(), gs.out.collect(), plate_ids )
	emit: 
		doc.out
}

workflow QC_AND_DOC_REPEAT {
	take: 
		est_counts_file 
		pc_pseudoalign_file
		plate_ids
	main:
		qc( est_counts_file, pc_pseudoalign_file )
	  gs(qc.out.collect(), plate_ids)
	  doc( qc.out.collect(), gs.out.collect(), plate_ids )
	emit: 
		doc.out
}

workflow QC_AND_DOC_plate_all {
	take: 
		est_counts_file 
		pc_pseudoalign_file
	main:
		def plate_ids = Channel.fromList( ['all'] )
		qc( est_counts_file, pc_pseudoalign_file )
	  gs(qc.out.collect(), plate_ids)
	  doc( qc.out.collect(), gs.out.collect(), plate_ids )
	emit: 
		doc.out
}

workflow {
  def read_pairs = Channel.fromFilePairs("${params.reads}/${params.pattern}.fastq.gz")
  
  quantification(read_pairs)  
  
  p_kal( quantification.out.quants.collect() )
	
  //def plate_ids = Channel.fromList( ['BC003BLCFUQ4','BC003BLHSCQ1','BC003BLMkP38Q2','BC003BLMkP38Q3'] )
  def plate_ids = Channel.fromList(params.plate_ids)

  //plate_ids.view()
  q_merge(p_kal.out, plate_ids )

	QC_AND_DOC(q_merge.out.estcounts, p_kal.out, plate_ids)
	//QC_AND_DOC.out.toList().view()
	
  mat_merge( QC_AND_DOC.out.toList() )
  
  tx2g( mat_merge.out )
	
	if(params.plate_ids.size() > 1) {
		QC_AND_DOC_plate_all(mat_merge.out, p_kal.out)
	}
	
	
}






