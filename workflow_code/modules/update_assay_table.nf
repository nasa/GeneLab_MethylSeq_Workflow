process UPDATE_ASSAY_TABLE {

    input:
        val(meta)
        path(runsheet)
        path(isa_zip)
        path(raw_multiqc_zip)
    
    output:
        path("a_*.txt"), emit: assay_table, optional: true

    script:
    def is_deduped = (meta.rrbs || params.skip_dedupe) ? "" : "--is_deduplicated"
    """
    update_assay_table.py \\
            --assay_suffix ${meta.assay_suffix} \\
            --runsheet ${runsheet} \\
            --glds_accession ${params.post_processing.glds_accession} \\
            --isa_zip ${isa_zip} \\
            ${is_deduped} \\
            --read_counts_from_multiqc ${raw_multiqc_zip}
    """
}