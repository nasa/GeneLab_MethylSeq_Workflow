process CLEAN_OUTPUT_PATHS { 
    tag "Purging paths of ${ file(alignment_report).name }"

    input:
        tuple path(alignment_report), path(outdir)

    output:
        path(alignment_report)


    script:
    """
    # Clean Alignment Reports in place
    clean_paths.sh "${alignment_report}"

    # Overwrite file in original folder
    cp "${alignment_report}" "${outdir}"
    """
}