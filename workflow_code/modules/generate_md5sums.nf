process GENERATE_MD5SUMS {
  // Generates tabular data for GeneLab raw (if needed) and processed data

    input:
        val(assay_suffix)
        val(processed_dir)
        val(done_token)
    
    output:
        path("raw_md5sum${assay_suffix}.tsv"), emit: raw_md5sum, optional: true
        path("processed_md5sum${assay_suffix}.tsv"), emit: processed_md5sum

    script:
        def raw_md5sum_flag = params.post_processing.include_raw_data ? "--generate_raw_md5sums" : ''
        """
        generate_md5sums.py --outdir ${processed_dir} --assay_suffix ${assay_suffix} ${raw_md5sum_flag}
        """
}