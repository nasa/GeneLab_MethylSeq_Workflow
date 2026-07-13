process TRIMGALORE {
  tag "${ meta.id }"

  input:
    tuple val(meta), path(reads)

  output:
    tuple val(meta), path("${ meta.id }*trimmed.fastq.gz"), emit: reads
    path("${ meta.id }*.txt"), emit: reports
    path("versions.yml"), emit: version

  script:
    def preset = params.library_presets[meta.kit]
    // Warn if no preset found and all clip params are null
    if (!preset) {
      def relevantParams = meta.paired_end ?
        [params.clip_r1, params.clip_r2, params.three_prime_clip_r1, params.three_prime_clip_r2] :
        [params.clip_r1, params.three_prime_clip_r1]
      if (relevantParams.every { it == null }) {
        log.warn "Library '${meta.kit}' not found in library_presets config and clip params were kept null, proceeding without any clipping.\n" +
            "If clipping is required, either pass clip params directly (--clip_r1, --three_prime_clip_r1${meta.paired_end ? ', --clip_r2, --three_prime_clip_r2' : ''}) or add a preset for '${meta.kit}' to library_presets in your config."
      }
      preset = [:]
    }

    // Helper closure to build optional clip flag, prioritizing user param over preset
    def clipFlag = { paramValue, presetValue, flagName ->
      if (paramValue != null) return "--${flagName} ${paramValue}"
      if (presetValue != null) return "--${flagName} ${presetValue}"
      return '' // nothing passed
    }

    // Resolve flags
    def clipR1Flag = clipFlag(params.clip_r1, preset.clip_r1, 'clip_R1')
    def clipR2Flag = meta.paired_end ? clipFlag(params.clip_r2, preset.clip_r2, 'clip_R2') : ''
    def tpcR1Flag  = clipFlag(params.three_prime_clip_r1, preset.three_prime_clip_r1, 'three_prime_clip_R1')
    def tpcR2Flag  = meta.paired_end ? clipFlag(params.three_prime_clip_r2, preset.three_prime_clip_r2, 'three_prime_clip_R2') : ''

    def assay_suffix = params.assay_suffix ? params.assay_suffix : "${meta.assay_suffix}"
    def rrbs_flag = (meta.rrbs && meta.kit != "nugen") ? '--rrbs' : ''
    def non_directional_flag = (rrbs_flag && (params.non_directional || preset.non_directional)) ? '--non_directional' : ''

    """
    trim_galore --gzip \
    --cores $task.cpus \
    --phred33 \
    --output_dir . \
    ${ meta.paired_end ? '--paired' : '' } \
    ${ meta.kit == "nugen" ? '-a AGATCGGAAGAGC' : '' } \
    ${ meta.kit == "nugen" && meta.paired_end ? '-a2 AAATCAAAAAAAC' : '' } \
    ${ rrbs_flag } \
    ${ non_directional_flag } \
    ${clipR1Flag} \
    ${clipR2Flag} \
    ${tpcR1Flag} \
    ${tpcR2Flag} \
    ${reads}

    # rename trimmed fastq files with assay suffix and _trimmed suffix
    ${ meta.paired_end ? \
      "mv ${ meta.id }${ assay_suffix }_R1_raw_val_1.fq.gz ${ meta.id }${ assay_suffix }_R1_trimmed.fastq.gz; \
      mv ${ meta.id }${ assay_suffix }_R2_raw_val_2.fq.gz ${ meta.id }${ assay_suffix }_R2_trimmed.fastq.gz" : \
      "mv ${ meta.id }${ assay_suffix }_raw_trimmed.fq.gz ${ meta.id }${ assay_suffix }_trimmed.fastq.gz" }

    # rename trimming report files
    ${ meta.paired_end ? \
      "mv ${ meta.id }${ assay_suffix }_R1_raw.fastq.gz_trimming_report.txt ${ meta.id }_R1${ assay_suffix }_trimming_report.txt; \
      mv ${ meta.id }${ assay_suffix }_R2_raw.fastq.gz_trimming_report.txt ${ meta.id }_R2${ assay_suffix }_trimming_report.txt" : \
      "mv ${ meta.id }${ assay_suffix }_raw.fastq.gz_trimming_report.txt ${ meta.id }${ assay_suffix }_trimming_report.txt"}

    echo '"${task.process}":' > versions.yml
    echo "    TrimGalore!: \$(trim_galore -v | sed -n 's/.*version //p' | head -n 1)" >> versions.yml
    echo "    Cutadapt: \$(cutadapt --version)" >> versions.yml
    """
}