process QUALIMAP {
  tag "${ meta.id }"

  input:
    path(genomeGtf)
    tuple val(meta), path(sorted_bam)

  output:
    tuple val(meta), path("${ meta.id }*_qualimap"), emit: qualimap_dir
    path  "versions.yml", emit: version

  script:
    def assay_suffix = params.assay_suffix ? params.assay_suffix : "${meta.assay_suffix}"
    def bismark_suffix = meta.rna ? "_hisat2" : "_bt2"
    """
    qualimap bamqc -bam ${ sorted_bam } \
        -gff ${genomeGtf} \
        -outdir ${ meta.id }${ assay_suffix }_bismark${ bismark_suffix }_qualimap/ \
        --collect-overlap-pairs \
        --java-mem-size=8G \
        -nt ${task.cpus}

    echo '"${task.process}":' > versions.yml
    echo "    qualimap: \$(qualimap 2>&1 | sed -n 's/^.*QualiMap v\\.\\([0-9.]\\+\\).*\$/\\1/p')" >> versions.yml
    """
}

process ZIP_QUALIMAP {
    tag "${ meta.id }"
    label "zip"
 
    input:
        tuple val(meta), path(qualimap_dir)

    output:
        path("${qualimap_dir}.zip"), emit: qualimap_zip
        path("versions.yml"), emit: version

    script:
        """
        # Zipping
        zip -q -r ${qualimap_dir}.zip ${qualimap_dir}
        zip -h | grep "Zip" | sed -E 's/(Zip.+\\)).+/\\1/' > versions.yml
        """
}
