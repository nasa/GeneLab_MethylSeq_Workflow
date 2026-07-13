process VALIDATE_PROCESSING {
    
    input:
        val meta
        path runsheet
        val processed_dir
        val done_token


    output:
        path("*methylseq-validation.log"), emit: validation_log

    script:
      def is_single = meta.paired_end ? "" : "--single_ended"
      def is_deduped = (meta.rrbs || params.skip_dedupe) ? "" : "--deduped"
      def raw_fastq = params.post_processing.include_raw_data ? "--include_raw_fastq" : ""
      """    
      GL-validate-processed-methylseq-data.py \\
            --assay_suffix '${meta.assay_suffix}' \\
            --runsheet '${runsheet}' \\
            --outdir '${processed_dir}' \\
            --GLDS_ID '${params.post_processing.glds_accession}' \\
            ${is_single} \\
            ${is_deduped} \\
            ${raw_fastq}
      """
}