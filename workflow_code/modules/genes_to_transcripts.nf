process GENES_TO_TRANSCRIPTS {
  storeDir "${ params.derived_store_path }/Genome_GTF_BED_Files/${ref_source}/${ref_source.toLowerCase().contains('ensembl') ? ref_version + '/' : ''}${organism_sci}"

  input:
    path(genome_gtf)
    val(organism_sci) // Used for defining storage location
    val(ref_source) // Used for defining storage location
    val(ref_version) // Used for defining storage location

  output:
    path("${ genome_gtf.baseName }-gene-to-transcript-map.tsv")

  script:
  """
  awk ' \$3 == \"transcript\" ' ${ genome_gtf } | cut -f 9 | tr -s \";\" \"\t\" | \
    cut -f 1,3 | tr -s \" \" \"\t\" | cut -f 2,4 | tr -d '\"' \
    > ${ genome_gtf.baseName }-gene-to-transcript-map.tsv
  """
}
