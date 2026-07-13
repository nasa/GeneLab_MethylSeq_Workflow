#!/usr/bin/env nextflow

// Declare syntax version
nextflow.enable.dsl=2

def colorCodes = [
    c_line: "┅" * 70,
    c_back_bright_red: "\u001b[41;1m",
    c_bright_green: "\u001b[32;1m",
    c_blue: "\033[0;34m",
    c_yellow: "\u001b[33;1m",
    c_reset: "\033[0m"
]

// making sure specific expected nextflow version is being used
if( ! nextflow.version.matches( workflow.manifest.nextflowVersion ) ) {
    println "\n    This workflow requires Nextflow version $workflow.manifest.nextflowVersion, but version $nextflow.version is currently active."
    println "\n    You can set the proper version for this terminal session by running: `export NXF_VER=$workflow.manifest.nextflowVersion`"

    println "\n  Exiting for now.\n"
    exit 1
}

////////////////////////////////////////////////////
/* --              DEBUG WARNING               -- */
////////////////////////////////////////////////////
if ( params.truncate_to ) {

    println( "\n    WARNING: DEBUG OPTIONS ARE ENABLED!\n" )

    params.truncate_to ? println( "        - Truncating reads to first ${ params.truncate_to } records" ) : null
    params.limit_samples_to ? println( "        - Limiting number of samples to process to ${ params.limit_samples_to }" ) : null
    params.force_single_end ? println( "        - Forcing single-end processing (using read 1 only)" ) : null

    println ""

}


////////////////////////////////////////////////////
/* --                WORKFLOW                  -- */
////////////////////////////////////////////////////

include { paramsHelp } from 'plugin/nf-schema'

include { STAGE_ANALYSIS } from './subworkflows/stage_analysis.nf'

include { FASTQC as RAW_FASTQC } from './modules/fastqc.nf'
include { FASTQC as TRIMMED_FASTQC } from './modules/fastqc.nf'
include { TRIMGALORE } from './modules/trimgalore.nf'
include { NUGEN } from './modules/nugen.nf'

include { PARSE_ANNOTATIONS_TABLE } from './modules/parse_annotations_table.nf'
include { DOWNLOAD_REFERENCES } from './modules/download_references.nf'

include { BUILD_BISMARK; ALIGN_BISMARK; DEDUPE; EXTRACT_CALLS; BISMARK_REPORT; BISMARK_SUMMARY } from './modules/bismark.nf'
include { QUALIMAP; ZIP_QUALIMAP } from './modules/qualimap.nf'
include { DIFF_METHYLATION } from './modules/methylation.nf'

include { SAMTOOLS_SORT as SORT_BAM } from './modules/samtools.nf'
include { SAMTOOLS_SORT as SORT_DEDUPED_BAM } from './modules/samtools.nf'

include { TO_PRED; TO_BED } from './modules/ucsc.nf'
include { GENES_TO_TRANSCRIPTS } from './modules/genes_to_transcripts.nf'

include { MULTIQC as RAW_MULTIQC } from './modules/multiqc.nf'
include { MULTIQC as TRIMMED_MULTIQC } from './modules/multiqc.nf'
include { MULTIQC as ALIGN_MULTIQC } from './modules/multiqc.nf'

include { SOFTWARE_VERSIONS } from './modules/software_versions.nf'
include { GENERATE_PROCESSING_PROTOCOL } from './modules/generate_protocols.nf'

ch_dp_tools_plugin = params.dp_tools_plugin ? channel.value(file(params.dp_tools_plugin)) : channel.value(file("$projectDir/bin/dp_tools__methylseq_dna"))

ch_runsheet = params.runsheet_path ? channel.fromPath(params.runsheet_path) : null
ch_isa_archive = params.isa_archive_path ? channel.fromPath(params.isa_archive_path) : null

ch_multiqc_config = params.multiqc_config ? channel.fromPath( params.multiqc_config ) : channel.fromPath("NO_FILE")

workflow {
    main:
        
    // Stage analysis setup (inputs, and raw reads)
    STAGE_ANALYSIS(
        ch_dp_tools_plugin,
        params.accession,
        ch_isa_archive,
        ch_runsheet,
        params.api_url
    )
    samples = STAGE_ANALYSIS.out.samples
    raw_reads = STAGE_ANALYSIS.out.raw_reads
    runsheet_path = STAGE_ANALYSIS.out.runsheet_path
    isa_archive = STAGE_ANALYSIS.out.isa_archive
    osd_accession = STAGE_ANALYSIS.out.osd_accession
    glds_accession = STAGE_ANALYSIS.out.glds_accession
    ch_dp_tools_version = STAGE_ANALYSIS.out.dp_tools_version

    // Extract sample IDs to samples.txt which is used in multiqc processes
    samples | map { it[0].id }
            | collectFile(name: "samples.txt", sort: true, newLine: true, storeDir: "processing_scripts")
            | set { ch_samples_txt }

    // Extract metadata from the first sample and set it as a channel
    samples | first
            | map { meta, reads -> meta }
            | set { ch_meta }

    // Extract assay suffix from metadata to be used in processes where meta is not passed
    ch_meta | map { meta -> meta.assay_suffix }
            | set { assay_suffix }
    
    RAW_FASTQC( raw_reads )
    RAW_FASTQC.out.fastqc | map { it -> [ it[1], it[2] ] } 
                          | flatten 
                          | collect 
                          | set { raw_mqc_ch }
    RAW_MULTIQC( ch_samples_txt, raw_mqc_ch, ch_multiqc_config, "raw", assay_suffix )

    TRIMGALORE( raw_reads )

    // Run NUGEN-specific trimming if kit is nugen, otherwise pass through trimmed reads from TRIMGALORE
    NUGEN( TRIMGALORE.out.reads )
    ch_trimmed_reads = NUGEN.out.reads.mix( 
                            TRIMGALORE.out.reads.filter { meta, reads -> meta.kit != "nugen" }
                        )
    
    TRIMMED_FASTQC( ch_trimmed_reads )
    TRIMMED_FASTQC.out.fastqc | map { it -> [ it[1], it[2] ] } 
                              | flatten 
                              | concat( TRIMGALORE.out.reports ) 
                              | collect 
                              | set { trimmed_mqc_ch }
    TRIMMED_MULTIQC( ch_samples_txt, trimmed_mqc_ch, ch_multiqc_config, "trimmed", assay_suffix )
    
    // Extract organism name from metadata
    ch_meta | map { meta -> meta.organism_sci }
        | set { organism_sci }
    
    // Use reference input and gene annotations file workflow params if provided
    if ( params.ref_fasta && params.ref_gtf ) {
        channel.value( params.ref_source ) | set { reference_source }
        channel.value( params.ref_version ) | set { reference_version }
        channel.value( params.ref_fasta ) | set { reference_fasta_url }
        channel.value( params.ref_gtf ) | set { reference_gtf_url }
        channel.value( params.gene_annotations_file ) | set { gene_annotations_url }
    } else{
        // Use annotations table to get reference inputs, organism-specific gene annotations file
        PARSE_ANNOTATIONS_TABLE( params.reference_table, organism_sci )
        reference_source = PARSE_ANNOTATIONS_TABLE.out.reference_source
        reference_version = PARSE_ANNOTATIONS_TABLE.out.reference_version
        reference_fasta_url = PARSE_ANNOTATIONS_TABLE.out.reference_fasta_url
        reference_gtf_url = PARSE_ANNOTATIONS_TABLE.out.reference_gtf_url
        gene_annotations_url = PARSE_ANNOTATIONS_TABLE.out.gene_annotations_url
     }

    DOWNLOAD_REFERENCES( organism_sci, reference_fasta_url, reference_gtf_url, reference_source, reference_version )
    genome_references = DOWNLOAD_REFERENCES.out.reference_files
    
    BUILD_BISMARK( 
      genome_references,
      reference_source,
      reference_version,
      assay_suffix,
      organism_sci
    )

    ALIGN_BISMARK( ch_trimmed_reads, BUILD_BISMARK.out.build )
    SORT_BAM( ALIGN_BISMARK.out.bam, ".bam" )

    QUALIMAP( genome_references | map { it[1] }, SORT_BAM.out.sorted_bam )
    ZIP_QUALIMAP( QUALIMAP.out.qualimap_dir )

    // Dedupe only if dataset is not RRBS
    DEDUPE ( ALIGN_BISMARK.out.bam )
    SORT_DEDUPED_BAM( DEDUPE.out.bam, ".deduplicated.bam" )

    // Mix both channels - DEDUPE runs when NOT rrbs AND NOT skip_dedupe
    // For RRBS samples or when skip_dedupe is true, use original BAM from ALIGN_BISMARK
    ch_aligned_reads = DEDUPE.out.bam.mix(
                            ALIGN_BISMARK.out.bam.filter { meta, bam -> meta.rrbs || params.skip_dedupe }
                        )    
    
    ch_dedupe_report = DEDUPE.out.report.mix(
                            ALIGN_BISMARK.out.bam
                                .filter { meta, bam -> meta.rrbs || params.skip_dedupe }
                                .map { meta, bam -> [ meta, file("NO_FILE") ] }
                        )

    EXTRACT_CALLS( ch_aligned_reads, BUILD_BISMARK.out.build )

    // Compiling all individual sample reports in one channel to send to BISMARK_REPORT
    ch_all_sample_reports = ALIGN_BISMARK.out.report | join( ALIGN_BISMARK.out.stats )
                                                     | join( EXTRACT_CALLS.out.splitting_report )
                                                     | join( EXTRACT_CALLS.out.bias )
                                                     | join( ch_dedupe_report )
                                                     | map { meta, align_report, align_stats, splitting_report, bias, dedupe_report ->
                                                            // Replace NO_FILE with empty list 
                                                            def dedupe = dedupe_report.name == "NO_FILE" ? [] : dedupe_report
                                                            [ meta, align_report, align_stats, splitting_report, bias, dedupe ]
                                                     }
    BISMARK_REPORT( ch_all_sample_reports )

    // Creating a channel holding all input files for bismark2summary (bam files, align reports, splitting reports, dedupe reports if they exist)
    ch_bams_and_all_reports = ALIGN_BISMARK.out.bam | join( ch_all_sample_reports ) 
                                                    | map { it -> it[ 1..it.size() - 1 ] } // Drop meta, keep all reports and bam in one tuple for BISMARK_SUMMARY
                                                    | collect
    BISMARK_SUMMARY( ch_bams_and_all_reports, assay_suffix )

    align_mqc_ch = QUALIMAP.out.qualimap_dir | join( ch_all_sample_reports ) 
                                             | map { it -> 
                                                    // Drop meta, then filter out empty lists in case of missing dedupe reports
                                                    it[1..-1].findAll { elem -> elem != [] && elem } 
                                             } 
                                             | collect
    ALIGN_MULTIQC( ch_samples_txt, align_mqc_ch, ch_multiqc_config, "align_and_bismark", assay_suffix  )

    TO_PRED( 
        genome_references | map { it[1] }, 
        organism_sci,
        reference_source,
        reference_version
    )

    TO_BED( 
        TO_PRED.out.genome_pred, 
        organism_sci,
        reference_source,
        reference_version
    )

    GENES_TO_TRANSCRIPTS( 
        genome_references | map { it[1] }, 
        organism_sci,
        reference_source,
        reference_version
    )

    DIFF_METHYLATION(
        EXTRACT_CALLS.out.cov | map { meta, cov -> cov } | collect,
        runsheet_path,
        organism_sci,
        TO_BED.out.genome_bed,
        GENES_TO_TRANSCRIPTS.out,
        gene_annotations_url,
        reference_source,
        reference_version,
        assay_suffix
    )
    
    // Software Version Capturing
        ch_software_versions = channel.empty()

        ch_dp_tools_version = ch_dp_tools_version ? ch_dp_tools_version : channel.empty()

        nf_version = '"NEXTFLOW":\n    nextflow: '.concat("${nextflow.version}\n")
        ch_nextflow_version = channel.value(nf_version)

        // Mix in versions from each process
        ch_software_versions = ch_software_versions
            | mix(ch_dp_tools_version)
            | mix(RAW_FASTQC.out.version)
            | mix(TRIMGALORE.out.version)
            | mix(RAW_MULTIQC.out.version)
            | mix(ALIGN_BISMARK.out.version)
            | mix(SORT_BAM.out.version)  
            | mix(QUALIMAP.out.version)
            | mix(TO_PRED.out.version)
            | mix(TO_BED.out.version)
            | mix(DIFF_METHYLATION.out.version)
            | mix(ch_nextflow_version)
        
        // Process the versions:
        ch_software_versions 
            | unique  
            | collectFile(
                newLine: true, 
                cache: false
            )
            | set { ch_final_software_versions }
        
        // Convert software versions combined yaml to markdown table
        SOFTWARE_VERSIONS( ch_final_software_versions, assay_suffix )

    // Generating protocol
    GENERATE_PROCESSING_PROTOCOL(
        ch_meta,
        SOFTWARE_VERSIONS.out.software_versions_yaml,
        genome_references,
        reference_source,
        reference_version
    )

    publish:
    // Metadata
    isa_archive = isa_archive
    runsheet = runsheet_path

    // Raw reads
    raw_reads = raw_reads

    // FastQC
    raw_fastqc = RAW_FASTQC.out.fastqc
    trimmed_fastqc = TRIMMED_FASTQC.out.fastqc

    // Trimmed reads and reports
    trimmed_reads = ch_trimmed_reads // TRIMGALORE.out.reads or NUGEN.out.reads depending on kit
    trimgalore_reports = TRIMGALORE.out.reports

    // MultiQC
    raw_multiqc_data = RAW_MULTIQC.out.data_folder
    raw_multiqc_html = RAW_MULTIQC.out.html
    trimmed_multiqc_data = TRIMMED_MULTIQC.out.data_folder
    trimmed_multiqc_html = TRIMMED_MULTIQC.out.html
    align_multiqc_data = ALIGN_MULTIQC.out.data_folder
    align_multiqc_html = ALIGN_MULTIQC.out.html

    // Alignment reports
    alignment_reports = ALIGN_BISMARK.out.report
    alignment_stats = ALIGN_BISMARK.out.stats
    sorted_bams = SORT_BAM.out.sorted_bam
    sorted_bais = SORT_BAM.out.sorted_bai
    dedupe_reports = ch_dedupe_report
                        .filter { meta, report -> report.name != "NO_FILE" } // Filter out samples where dedupe was not run and report is just a placeholder
    deduped_sorted_bams = SORT_DEDUPED_BAM.out.sorted_bam
    deduped_sorted_bais = SORT_DEDUPED_BAM.out.sorted_bai
    qualimap = ZIP_QUALIMAP.out.qualimap_zip

    // Methylation calls and reports
    methylation_contexts = EXTRACT_CALLS.out.contexts
    methylation_bedgraph = EXTRACT_CALLS.out.bedGraph
    methylation_cov = EXTRACT_CALLS.out.cov
    methylation_cytosine_report = EXTRACT_CALLS.out.cytosine_report
    methylation_bias = EXTRACT_CALLS.out.bias
    methylation_splitting_report = EXTRACT_CALLS.out.splitting_report
    methylation_cytosine_summary = EXTRACT_CALLS.out.cytosine_summary
    bismark_report = BISMARK_REPORT.out.html
    bismark_summary_txt = BISMARK_SUMMARY.out.summary_txt
    bismark_summary_html = BISMARK_SUMMARY.out.summary_html

    // Differential methylation analysis
    sample_table = DIFF_METHYLATION.out.sample_table
    contrasts = DIFF_METHYLATION.out.contrasts
    diff_methylation_bases = DIFF_METHYLATION.out.bases
    diff_methylation_tiles = DIFF_METHYLATION.out.tiles

    // GeneLab
    software_versions = SOFTWARE_VERSIONS.out.software_versions_md
    processing_protocol = GENERATE_PROCESSING_PROTOCOL.out.proc_protocol

}

output {
    // Metadata
    isa_archive {
        path "Metadata"
    }

    runsheet {
        path "Metadata"
    }

    // Raw reads
    raw_reads {
        path { meta, reads -> "Merged_Sequence_Data" }
    }

    // FastQC
    raw_fastqc {
        path { meta, html, zip -> "Merged_Sequence_Data/FastQC_Reports" }
    }

    trimmed_fastqc {
        path { meta, html, zip -> "Trimmed_Sequence_Data/FastQC_Reports" }
    }

    // Trimmed reads and reports
    trimmed_reads {
        path { meta, reads -> "Trimmed_Sequence_Data" }
    }

    trimgalore_reports {
        path "Trimmed_Sequence_Data/Trimming_Reports"
    }

    // MultiQC
    raw_multiqc_data {
        path "Merged_Sequence_Data/MultiQC_Reports"
    }

    raw_multiqc_html {
        path "Merged_Sequence_Data/MultiQC_Reports"
    }

    trimmed_multiqc_data {
        path "Trimmed_Sequence_Data/MultiQC_Reports"
    }

    trimmed_multiqc_html {
        path "Trimmed_Sequence_Data/MultiQC_Reports"
    }

    align_multiqc_data {
        path "Methylation_Call_Data/MultiQC_Reports"
    }

    align_multiqc_html {
        path "Methylation_Call_Data/MultiQC_Reports"
    }

    // Alignment
    alignment_reports {
        path { meta, report -> "Aligned_Sequence_Data/Alignment_Reports" }
    }

    alignment_stats {
        path { meta, stats -> "Aligned_Sequence_Data/Alignment_Reports" }
    }

    sorted_bams {
        path { meta, bam -> "Aligned_Sequence_Data" }
    }

    sorted_bais {
        path "Aligned_Sequence_Data"
    }

    dedupe_reports {
        path { meta, report -> "Aligned_Sequence_Data/Deduplication_Reports" }
    }

    deduped_sorted_bams {
        path { meta, bam -> "Aligned_Sequence_Data" }
    }

    deduped_sorted_bais {
        path "Aligned_Sequence_Data"
    }

    qualimap {
        path "Aligned_Sequence_Data/Alignment_Reports"
    }

    // Methylation calls
    methylation_contexts {
        path { meta, files -> "Methylation_Call_Data" }
    }

    methylation_bedgraph {
        path { meta, file -> "Methylation_Call_Data" }
    }

    methylation_cov {
        path { meta, file -> "Methylation_Call_Data" }
    }

    methylation_cytosine_report {
        path { meta, file -> "Methylation_Call_Data" }
    }

    methylation_bias {
        path { meta, file -> "Methylation_Call_Data/Methylation_Call_Reports" }
    }

    methylation_splitting_report {
        path { meta, file -> "Methylation_Call_Data/Methylation_Call_Reports" }
    }

    methylation_cytosine_summary {
        path { meta, file -> "Methylation_Call_Data/Methylation_Call_Reports" }
    }

    bismark_report {
        path { meta, html -> "Methylation_Call_Data/Methylation_Call_Reports" }
    }

    bismark_summary_txt {
        path "Methylation_Call_Data/Methylation_Call_Reports/Methylation_Call_Summary_Reports"
    }

    bismark_summary_html {
        path "Methylation_Call_Data/Methylation_Call_Reports/Methylation_Call_Summary_Reports"
    }

    // Differential methylation
    sample_table {
        path "Differential_Methylation_Analysis_Data"
    }

    contrasts {
        path "Differential_Methylation_Analysis_Data"
    }

    diff_methylation_bases {
        path "Differential_Methylation_Analysis_Data"
    }

    diff_methylation_tiles {
        path "Differential_Methylation_Analysis_Data"
    }

    // GeneLab
    software_versions {
        path "GeneLab"
    }

    processing_protocol {
        path "GeneLab"
    }
}
