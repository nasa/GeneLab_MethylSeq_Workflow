process CLEAN_MULTIQC_PATHS { 
    tag "Purging paths of ${ file(multiqc_data).name }"

    input:
        tuple path(multiqc_data), path(outdir)

    output:
        path("${outdir}/*.zip")

    script:
    """
    clean_multiqc_paths.py "${multiqc_data}" "${outdir}"
    """
}