// standard run case:'nextflow run scqc_nf.sh -c PIP-3144.config -with-dag flowchart.png -with-report -resume'

nextflow.enable.dsl=2

process quantification {
	debug true
	label "image_rscater"
  publishDir "$params.quantificationsoutdir"
  
	
	tag "$sampleId"
  input:
  	tuple val(sampleId), path(R1), path(R2)

  output:
  	path "${sampleId}" , emit: quants

  """
  kallisto quant -t 1 -i ${params.idx} -o $sampleId -b 100 $R1 $R2
  """
	
  /*
	tag "${sample_id}"
	input:
		val sample_id 

	output: 
		path "${sample_id}", emit: quants

	"""
  kallisto quant -t 2 -i ${params.idx} -o $sample_id -b 100 $sample_id
  """
  */
}

process p_kal {
	label "image_rscater"

	// output percent_pseudoaligned.txt to qc_dir/ by copying
	publishDir "$params.qcoutdir", mode: 'copy'
	
	input:
		// wait for all kallisto processes to finish 
		path quants_check 

	output:
		path "percent_pseudoaligned.txt"

	"""
	Rscript ${params.RScript_dir}kallisto_mapping_scrape.R ${params.quantificationsoutdir} percent_pseudoaligned.txt p_pseudoaligned
	"""
}


process q_merge {
	label "image_rscater"
	tag "$plate_id"	
	errorStrategy 'finish'
	publishDir "$params.qcoutdir/$plate_id", mode: 'copy'

	input:
		// wait for p_kal to output percent_pseudoaligned.txt - will be ln to current process dir, although it is not used.
		path p_kal
		val plate_id 
	
	output:
		path "est_counts${plate_id}_matrix.tsv" , emit: estcounts
		val plate_id, emit: plate_id

	// R script generate, e.g., $workdir/est_countsCU5DAY0_matrix.tsv
	// and copy to qc_dir/$plate_id location
	"""	
	workdir="\$(pwd)" && \
	Rscript ${params.RScript_dir}merge_kallisto_quant.R ${params.quantificationsoutdir} est_counts ${plate_id} \${workdir};
	"""
}

process qc {
	label "image_rscater"
	errorStrategy 'finish'
	beforeScript 'export HDF5_DISABLE_VERSION_CHECK=1'
	tag "$plate_id"
	
	// publishDir "$params.qcoutdir/$plate_id", mode: 'copy', pattern: "{*.tsv, *.pdf, *.Rdata}"
	publishDir "$params.qcoutdir/$plate_id", mode: 'copy'

	input:
		// e.g., quants_dir/est_countsCU5DAY0_matrix.tsv
		path est_counts_file 
		path pc_pseudoalign_file
		val plate_id

	output:
		path "qc_${est_counts_file}.done", emit: complete_check
		path "${plate_id}qc_for_doc.Rdata", emit: rdata
		path "QC_meanexp_vs_freq${plate_id}.pdf", emit: figpdf
		path "*.pdf"
		path "*.tsv"

	""" 
	workdir="\$(pwd)" && cd \${workdir} && \
	Rscript ${params.RScript_dir}scqc_from_matrix.meta.R ${est_counts_file} ${pc_pseudoalign_file} \${workdir} ${params.samplesheet} ${params.mtnamefile} && \
	touch qc_${est_counts_file}.done
	"""
}

process gs {
	errorStrategy 'finish'
	publishDir "$params.qcoutdir/$plate_id", mode: 'copy'
	label "image_rscater" 
	tag "$plate_id"

	input:
		val plate_id
		path figpdf

	output:
		path "rename_QC_meanexp_vs_freq${plate_id}.png"

	"""	
	gs -dNOPAUSE -dQUIET -dBATCH -sDEVICE=png16m -sOutputFile=rename_QC_meanexp_vs_freq${plate_id}.png -r256 \
    ${figpdf}
  """
}

process doc {
	label "image_rknit"
	tag "$plate_id"
	//maxForks 1

	publishDir "$params.qcoutdir/$plate_id", mode: 'copy'
	//cache false

	input:
		path qc_rdata
		path qc_complete_check
		path gs_png_file 
		val plate_id 
		path percent_pseudoaligned
	
	output:
		path "Finished_${plate_id}.txt", emit: complete_check
		path "${plate_id}_QC_report.pdf", emit: pdf

	"""
	workdir="\$(pwd)" && cd \${workdir} && echo \${workdir} && deref_rdata="\$(readlink ${qc_rdata})" && deref_png="\$(readlink ${gs_png_file})" && \
	Rscript -e \"options(warn=-1);objects<-\'\${deref_rdata}\';pngfile<-\'\${deref_png}\';mapping_file <- read.table(\'${percent_pseudoaligned}\');rmarkdown::render(\'${params.RScript_dir}QCreport.Rmd\', 'pdf_document', output_file=\'${plate_id}_QC_report.pdf\', output_dir=\'\${workdir}\')\" && \
	touch Finished_${plate_id}.txt
	"""
}

process mat_merge {
	// merge all plates est_counts${plate_id}_matrix.tsv into one
	publishDir "$params.qcoutdir/all", mode: 'copy'
  //cache false
  label "image_rscater"

	input:
		val est_counts_file_list 
 
	output:
		path 'est_countsall_plates.tsv'
	
	
	// singularity exec R_verse.v5.img Rscript plate_merge.R quants_dir/ \
	//	'Finished_CU5DAY0.txt Finished_CU7DAY0.txt Finished_CU5DAY7.txt Finished_CU7DAY7.txt';
	"""
	workdir="\$(pwd)" && cd \${workdir} && \
	Rscript ${params.RScript_dir}plate_merge.R \${workdir} est_counts \'${est_counts_file_list}\' nextflow;
	"""
}

process tx2g {
	publishDir "$params.qcoutdir/all", mode: 'copy'
	label "image_rscater"
	input:
		path all_plates_tsv 

	output:
		// path 'tx2g_finished.txt'
		path 'plates_as_genelevel.tsv'

	"""
	workdir="\$(pwd)" && cd \${workdir} && \
	Rscript ${params.RScript_dir}est_counts_tx2gene.R ${all_plates_tsv} \${workdir}/plates_as_genelevel.tsv ${params.species} ${params.trans2gen_tsv} 
	"""

}

workflow QC_AND_DOC {
	take: 
		est_counts_file 
		pc_pseudoalign_file
		plate_ids
	main:
		qc( est_counts_file, pc_pseudoalign_file, plate_ids )
	  gs( plate_ids, qc.out.figpdf )
	  doc( qc.out.rdata, qc.out.complete_check, gs.out, plate_ids, pc_pseudoalign_file )
	emit: 
		doc.out.complete_check
}

workflow QC_AND_DOC_plate_all {
	take: 
		est_counts_file 
		pc_pseudoalign_file
	main:
		def plate_ids = Channel.fromList( ['all'] )
		qc( est_counts_file, pc_pseudoalign_file, plate_ids )
		gs( plate_ids, qc.out.figpdf )
	  doc( qc.out.rdata, qc.out.complete_check, gs.out, plate_ids, pc_pseudoalign_file )
	emit: 
		doc.out.complete_check
}

workflow {
  def sample_ids = Channel.fromPath(params.samplesheet) | splitCsv(header:true) | map { row-> row.Sample_ID}
  
  def plate_ids = Channel.fromList(params.plate_ids)

	Channel
		.fromPath(params.samplecsv)
		.splitCsv(header:true)
		.map { row -> tuple(row.sampleId, file(row.R1), file(row.R2)) }
		.set {sample_ids_ch}

	quantification(sample_ids_ch)
  
  p_kal( quantification.out.quants.collect() )
	
  q_merge(p_kal.out, plate_ids )

	QC_AND_DOC(q_merge.out.estcounts, p_kal.out, q_merge.out.plate_id)
	
  mat_merge( q_merge.out.estcounts.collect() )
  
  tx2g( mat_merge.out )
	
	if(params.plate_ids.size() > 1) {
		QC_AND_DOC_plate_all(mat_merge.out, p_kal.out)
	}

	
}






