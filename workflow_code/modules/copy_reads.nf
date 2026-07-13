process COPY_READS {
    tag "${ meta.id }"

    input:
        tuple val(meta), path("?.gz")

    output:
        tuple val(meta), path("${meta.id}*.gz"), emit: raw_reads

    script:
        def assay_suffix = params.assay_suffix ? params.assay_suffix : "${meta.assay_suffix}"
        if ( meta.paired_end ) {
        """
        cp -P 1.gz ${meta.id}${assay_suffix}_R1_raw.fastq.gz
        cp -P 2.gz ${meta.id}${assay_suffix}_R2_raw.fastq.gz
        """
        } else {
        """
        cp -P 1.gz ${meta.id}${assay_suffix}_raw.fastq.gz
        """
        }
}
