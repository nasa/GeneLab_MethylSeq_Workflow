process MULTIQC {
    tag("Dataset-wide")
    
    input:
      path(sample_names)
      path("mqc_in/*") // any number of multiqc compatible files
      path(multiqc_config)
      val(mqc_label)
      val(suffix)
    
    output:
      path("${ mqc_label }_multiqc*_data"), emit: data_folder
      path("${ mqc_label }_multiqc*.html"), emit: html
      path("versions.yml"), emit: version

    script:
    def config_arg = multiqc_config.name != "NO_FILE" ? "--config ${multiqc_config}" : ""
    def assay_suffix = params.assay_suffix ? params.assay_suffix : "${suffix}"
    """
    multiqc \\
            --force \\
            --interactive \\
            -o . \\
            -n ${ mqc_label }_multiqc${ assay_suffix } \\
            ${ config_arg } \\
            .

    echo "${task.process}:" > versions.yml
    echo "    multiqc: \$(multiqc --version | sed -e "s/multiqc, version //g")" >> versions.yml
    """
}
