process GENERATE_README {
    
    input:
        val meta

    output:
        path("README${meta.assay_suffix}.txt"), emit: readme

    script:
      def is_paired = meta.paired_end ? '--paired-end' : ''
      def dedupe_skipped = (meta.rrbs || params.skip_dedupe) ? '--dedupe-skipped' : ''
      def raw_reads = params.post_processing.include_raw_data ? '--include-raw-reads' : ''
      """    
      GL-gen-processed-data-methylseq-readme.py \\
            --assay_suffix '${meta.assay_suffix}' \\
            ${is_paired} \\
            ${dedupe_skipped} \\
            --name '${params.post_processing.name}' \\
            --email '${params.post_processing.email}' \\
            --protocol-ID '${params.post_processing.protocol_id}' \\
            --osd-id '${params.post_processing.osd_accession}' \\
            ${raw_reads}
      """
}