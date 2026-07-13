process BUILD_BISMARK {
  tag "Refs: ${ genome_fasta }, ${ genome_gtf }, Source: ${ref_source}${ref_source.toLowerCase().contains('ensembl') ? ', Version: ' + ref_version : ''}"
  storeDir "${ params.derived_store_path }/BismarkIndices${suffix.toLowerCase().contains('rna') ? '_HT' : '_BT'}/${ref_source}/${ref_source.toLowerCase().contains('ensembl') ? ref_version + '/' : ''}${organism_sci}"

  label 'maxCPU'
  label 'big_mem'

  input:
    tuple path(genome_fasta), path(genome_gtf)
    val(ref_source)
    val(ref_version)
    val(suffix)
    val(organism_sci)

  output:
    val(bismark_index_dir), emit: build

  script:
    def aligner = suffix == "_GLRNAMethylSeq" ? "--hisat2" : "--bowtie2"
    bismark_index_dir = "${ params.derived_store_path }/BismarkIndices${suffix.toLowerCase().contains('rna') ? '_HT' : '_BT'}/${ref_source}/${ref_source.toLowerCase().contains('ensembl') ? ref_version + '/' : ''}${organism_sci}"
    """
    mkdir -p ${ bismark_index_dir }

    # Copy ref fasta to bismark index dir since it is required to be part of it by bismark
    cp ${ genome_fasta } ${ bismark_index_dir }

    bismark_genome_preparation ${aligner} --parallel ${task.cpus} ${ bismark_index_dir }
    bam2nuc --genome_folder ${ bismark_index_dir } --genomic_composition_only
    """
}

process ALIGN_BISMARK {
  tag "${ meta.id }"

  input:
    tuple val( meta ), path( reads )
    path(bismark_index_dir)

  output:
    tuple val(meta), path("${ meta.id }*.bam"), emit: bam
    tuple val(meta), path("${ meta.id }*nucleotide_stats.txt"), emit: stats
    tuple val(meta), path("${ meta.id }*report.txt"), emit: report
    path("versions.yml"), emit: version

  script:
    def aligner = meta.rna ? "--hisat2" : "--bowtie2"
    def bismark_suffix = meta.rna ? "_hisat2" : "_bt2"
    def preset = params.library_presets[meta.kit]
    if (!preset) {
      preset = [:]
    }
    def non_directional_flag = (params.non_directional || preset.non_directional) ? '--non_directional' : ''
    def pbat_flag = meta.kit == "pbat" ? '--pbat' : ''
    def assay_suffix = params.assay_suffix ? params.assay_suffix : "${meta.assay_suffix}"
    def input = meta.paired_end ? "-1 ${ reads[0] } -2 ${ reads[1] }" : "${ reads }"
    def aligner_version = meta.rna
      ? 'echo "    Hisat2: \$(hisat2 --version | awk \'{print \$3}\')" >> versions.yml'
      : 'echo "    Bowtie2: \$(bowtie2 --version | head -n1 | awk \'{print \$3}\')" >> versions.yml'
    """
    bismark ${aligner} \
      --bam \
      --parallel 4 \
      --gzip \
      --non_bs_mm \
      ${ non_directional_flag } \
      ${ pbat_flag } \
      --nucleotide_coverage \
      --genome_folder ${ bismark_index_dir } \
      ${ input }

    # remove R1/R2 and _trimmed in filename and add bismark suffix (_bt2 for DNA Methylation or _hisat2 for RNA Methylation)
    ${ meta.paired_end ? \
      "mv ${ meta.id }*_PE_report.txt ${ meta.id }${assay_suffix}_bismark${ bismark_suffix }_PE_report.txt; \
      mv ${ meta.id }*.bam ${ meta.id }${assay_suffix}_bismark${ bismark_suffix }_pe.bam" : \
      "mv ${ meta.id }*_SE_report.txt ${ meta.id }${assay_suffix}_bismark${ bismark_suffix }_SE_report.txt; \
      mv ${ meta.id }*.bam ${ meta.id }${assay_suffix}_bismark${ bismark_suffix }.bam" }

    mv ${ meta.id }*.nucleotide_stats.txt ${ meta.id }${assay_suffix}_bismark${ bismark_suffix }.nucleotide_stats.txt

    
    echo '"${task.process}":' > versions.yml
    echo "    Bismark: \$(bismark --version | sed -n 's/.*v//p')" >> versions.yml
    ${ aligner_version }
    """
}

process DEDUPE {
  tag "${ meta.id }"

  input:
    tuple val(meta), path(bam)

  output:
    tuple val(meta), path("${ meta.id }*deduplicated.bam"), emit: bam
    tuple val(meta), path("${ meta.id }*deduplication_report.txt"), emit: report

  when:
    !meta.rrbs && !params.skip_dedupe

  script:
    """
    deduplicate_bismark ${ bam }
    """
}

process EXTRACT_CALLS {
  tag "${ meta.id }"

  input:
    tuple val(meta), path(bam)
    path(bismark_index_dir)

  output:
    tuple val(meta), path("*_context_${ meta.id }*.txt.gz"), emit: contexts
    tuple val(meta), path("${ meta.id }*bedGraph.gz"), emit: bedGraph
    tuple val(meta), path("${ meta.id }*bismark.cov.gz"), emit: cov
    tuple val(meta), path("${ meta.id }*CpG_report.txt.gz"), emit: cytosine_report
    tuple val(meta), path("${ meta.id }*M-bias.txt"), emit: bias
    tuple val(meta), path("${ meta.id }*splitting_report.txt"), emit: splitting_report
    tuple val(meta), path("${ meta.id }*cytosine_context_summary.txt"), emit: cytosine_summary

  script:
    def preset = params.library_presets[meta.kit]
    if (!preset) {
      preset = [:]
    }

    // Helper closure to build optional ignore flag, prioritizing user param over preset
    def ignoreFlag = { paramValue, presetValue, flagName ->
      if (paramValue != null) return "--${flagName} ${paramValue}"
      if (presetValue != null) return "--${flagName} ${presetValue}"
      return '' // nothing passed
    }

    // Resolve flags
    def ignoreR1Flag = ignoreFlag(params.ignore_r1, preset.ignore_r1, 'ignore')
    def ignoreR2Flag = meta.paired_end ? ignoreFlag(params.ignore_r2, preset.ignore_r2, 'ignore_r2') : ''
    def ignore3primeR1Flag  = ignoreFlag(params.ignore_3prime_r1, preset.ignore_3prime_r1, 'ignore_3prime')
    def ignore3primeR2Flag  = meta.paired_end ? ignoreFlag(params.ignore_3prime_r2, preset.ignore_3prime_r2, 'ignore_3prime_r2') : ''

    // https://www.biostars.org/p/9594498/ - --genome_folder needs to be absolute path
    """
    bismark_methylation_extractor --parallel ${task.cpus} \
      --bedGraph \
      --gzip \
      --comprehensive \
      --cytosine_report \
      --genome_folder \${PWD}/${bismark_index_dir} \
      ${ ignoreR1Flag } \
      ${ ignore3primeR1Flag } \
      ${ ignoreR2Flag } \
      ${ ignore3primeR2Flag } \
      ${ bam }
    """
}

process BISMARK_REPORT {
  tag "${ meta.id }"

  input:
    tuple val(meta), path(align_report), path(align_stats), path(call_report), path(call_bias), path(dedupe_report)

  output:
    tuple val(meta), path("*_report.html"), emit: html

  script:
    def dedupe_flag = dedupe_report ? "--dedup_report ${dedupe_report}" : ""
    """
    bismark2report \
      --alignment_report ${ align_report } \
      --nucleotide_report ${ align_stats } \
      ${ dedupe_flag } \
      --splitting_report ${ call_report } \
      --mbias_report ${ call_bias }
    """
}

process BISMARK_SUMMARY {
  tag "Dataset-wide"

  input:
    file(unsorted_bams_and_reports)
    val(suffix)

  output:
    path("bismark_summary_report*.txt"), emit: summary_txt
    path("bismark_summary_report*.html"), emit: summary_html

  script:
    def assay_suffix = params.assay_suffix ? params.assay_suffix : "${suffix}"
    """
    bismark2summary -o "bismark_summary_report${ assay_suffix }" *bam
    """
}
