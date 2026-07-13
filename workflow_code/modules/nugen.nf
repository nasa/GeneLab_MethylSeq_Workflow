process NUGEN {
  // Stages the raw reads into appropriate publish directory
  publishDir "${params.outdir}/Trimmed_Sequence_Data",
      pattern: "*trimmed.fastq.gz",
      mode: params.publish_dir_mode
  tag "${ meta.id }"
  label 'low_cpu_med_memory'

  input:
    tuple val(meta), path("input/*")

  output:
    tuple val(meta), path("${ meta.id }*trimmed.fastq.gz"), emit: reads

  when:
    meta.kit == "nugen"

  script:
    def assay_suffix = params.assay_suffix ? params.assay_suffix : "${meta.assay_suffix}"
    if (meta.paired_end) {
    """
    trimRRBSdiversityAdaptCustomers.py -1 input/${ meta.id }${ assay_suffix }_R1_trimmed.fastq.gz -2 input/${ meta.id }${ assay_suffix }_R2_trimmed.fastq.gz

    mv input/${ meta.id }${ assay_suffix }_R1_trimmed.fastq_trimmed.fq.gz ${ meta.id }${ assay_suffix }_R1_trimmed.fastq.gz
    mv input/${ meta.id }${ assay_suffix }_R2_trimmed.fastq_trimmed.fq.gz ${ meta.id }${ assay_suffix }_R2_trimmed.fastq.gz
    """
    } else {
    """
    trimRRBSdiversityAdaptCustomers.py -1 input/${ meta.id }${ assay_suffix }_trimmed.fastq.gz

    mv input/${ meta.id }${ assay_suffix }_trimmed.fastq_trimmed.fq.gz ${ meta.id }${ assay_suffix }_trimmed.fastq.gz
    """
    }
    
}
