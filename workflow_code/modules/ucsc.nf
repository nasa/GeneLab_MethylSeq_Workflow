process TO_PRED {
  // Converts reference gtf into pred 
  storeDir "${ params.derived_store_path }/Genome_GTF_BED_Files/${ref_source}/${ref_source.toLowerCase().contains('ensembl') ? ref_version + '/' : ''}${organism_sci}"

  input:
    path(genome_gtf)
    val(organism_sci) // Used for defining storage location
    val(ref_source) // Used for defining storage location
    val(ref_version) // Used for defining storage location

  output:
    path("${ genome_gtf.baseName }.genepred"), emit: genome_pred
    path("versions.yml"), emit: version

  script: // https://github.com/nextflow-io/nextflow/issues/1359
  """
  gtfToGenePred ${ genome_gtf } ${ genome_gtf.baseName }.genepred

  echo '"${task.process}":' > versions.yml
  echo "    gtfToGenePred: 469" >> versions.yml
  """
}

process TO_BED {
  // Converts reference genePred into Bed format
  storeDir "${ params.derived_store_path }/Genome_GTF_BED_Files/${ref_source}/${ref_source.toLowerCase().contains('ensembl') ? ref_version + '/' : ''}${organism_sci}"

  input:
    path(genome_pred)
    val(organism_sci) // Used for defining storage location
    val(ref_source) // Used for defining storage location
    val(ref_version) // Used for defining storage location 

  output:
    path("${ genome_pred.baseName }.bed"), emit: genome_bed
    path("versions.yml"), emit: version

  script:
  """
  genePredToBed ${ genome_pred } ${ genome_pred.baseName }.bed

  echo '"${task.process}":' > versions.yml
  echo "    genePredToBed: 469" >> versions.yml
  """
}
