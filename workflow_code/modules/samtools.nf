process SAMTOOLS_SORT {
  tag "${ meta.id }"

  input:
    tuple val(meta), path(bam)
    val(ext)

  output:
    tuple val(meta), path("*_sorted${ ext }"), emit: sorted_bam
    path("*_sorted${ ext }.bai"), emit: sorted_bai
    path("versions.yml"), emit: version

  script:
    def assay_suffix = params.assay_suffix ? params.assay_suffix : "${meta.assay_suffix}"
    def bismark_suffix = meta.rna ? "_hisat2" : "_bt2"
    """
    samtools sort -@${task.cpus} \
        --write-index \
        --no-PG \
        -o ${ meta.id }${assay_suffix}_bismark${ bismark_suffix }_sorted${ ext }##idx##${ meta.id }${assay_suffix}_bismark${ bismark_suffix }_sorted${ ext }.bai \
        ${ bam }

    echo '"${task.process}":' > versions.yml
    echo "    samtools: \$(samtools --version | head -n1 | awk '{print \$2}')" >> versions.yml
    """
}
