process GENERATE_PROCESSING_PROTOCOL {

    input:
        val(ch_meta)
        path(software_versions_yaml)
        tuple path(reference_fasta), path(reference_gtf)
        val(ref_source)
        val(ref_version)

    output:
        path("processed_data_protocol*.txt"), emit: proc_protocol

    script:
        def mode = ch_meta.rna ? '--mode RNA' : ''
        def dedupe_skipped = (ch_meta.rrbs || params.skip_dedupe) ? '--dedupe_skipped' : ''
        def reference_source = ref_source ? "--reference_source ${ref_source}" : ''
        def reference_version = ref_version ? "--reference_version ${ref_version}" : ''
        def assay_suffix = params.assay_suffix ? params.assay_suffix : "${ch_meta.assay_suffix}"

        """
        generate_processed_protocol.py \
        ${mode} \
        --outdir . \
        --software_table ${software_versions_yaml} \
        --assay_suffix ${assay_suffix} \
        --workflow_version ${workflow.manifest.version} \
        --organism "${ch_meta.organism_sci}" \
        ${dedupe_skipped} \
        ${reference_source} \
        ${reference_version} \
        --reference_fasta ${reference_fasta} \
        --reference_gtf ${reference_gtf}
        """
}

process GENERATE_RAW_PROTOCOL {

    input:
        val(meta)
        path(software_versions_md)

    output:
        path("raw_data_protocol*.txt"), emit: raw_protocol

    when:
        params.post_processing.include_raw_data

    script:
        def raw_merged_flag = params.post_processing.merged_raw_data ? '--raw_merged' : ''

        """
        generate_raw_protocol.py \
        --outdir . \
        --software_table ${software_versions_md} \
        --assay_suffix ${meta.assay_suffix} \
        ${raw_merged_flag}
        """
}
