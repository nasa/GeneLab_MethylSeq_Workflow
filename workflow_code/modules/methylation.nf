process DIFF_METHYLATION {

  input:
    path("bismark_in/*")
    path(runsheet)
    val(organism_name)
    path(ref_bed)
    path(gene_transcript_map)
    val(annotation_file_link)
    val(ref_source)
    val(ref_version)
    val(suffix)

  output:
    path("SampleTable*.csv"), emit: sample_table
    path("contrasts*.csv"), emit: contrasts
    path("differential_methylation_bases*.csv"), emit: bases
    path("differential_methylation_tiles*.csv"), emit: tiles
    path("versions.yml"), emit: version

  script:
    def primary_keytype = ref_source == 'ensembl_plants' ? 'TAIR' : 'ENSEMBL'
    def assay_suffix = params.assay_suffix ? params.assay_suffix : "${suffix}"
    """
    differential_methylation.R \
        --bismark_methylation_calls_dir bismark_in \
        --path_to_runsheet ${runsheet} \
        --org_name ${organism_name.capitalize()} \
        --ref_bed_path ${ref_bed} \
        --gene_transcript_map_path ${gene_transcript_map}\
        --methylkit_output_dir . \
        --ref_ensemblVersion ${ref_version} \
        --ref_annotations_tab_link ${annotation_file_link} \
        --methRead_mincov ${params.methRead_mincov} \
        --getMethylDiff_difference ${params.MethylDiff_difference} \
        --getMethylDiff_qvalue ${params.MethylDiff_qvalue} \
        --tileMethylCounts_mincov ${params.tileMethylCounts_mincov} \
        --tileMethylCounts_winsize ${params.tileMethylCounts_winsize} \
        --tileMethylCounts_stepsize ${params.tileMethylCounts_stepsize} \
        --primary_keytype ${primary_keytype} \
        --file_suffix ${assay_suffix}

    echo '"${task.process}":' > versions.yml
    echo "    R: \$(Rscript -e 'cat(gsub(\" .*\",\"\",gsub(\"R version \",\"\",R.version\$version.string)))')" >> versions.yml
    echo "    Bioconductor: \$(Rscript -e \"cat(as.character(BiocManager::version()))\")" >> versions.yml
    echo "    methylKit: \$(Rscript -e \"cat(as.character(packageVersion('methylKit')))\" )" >> versions.yml
    echo "    genomation: \$(Rscript -e \"cat(as.character(packageVersion('genomation')))\" )" >> versions.yml
    echo "    tidyverse: \$(Rscript -e \"cat(as.character(packageVersion('tidyverse')))\" )" >> versions.yml
    echo "    dplyr: \$(Rscript -e \"cat(as.character(packageVersion('dplyr')))\" )" >> versions.yml
    """
}
