#!/usr/bin/env nextflow

// Declare syntax version
nextflow.enable.dsl=2

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

include { PARSE_RUNSHEET } from './subworkflows/parse_runsheet.nf'
include { CLEAN_MULTIQC_PATHS } from './modules/clean_multiqc_paths.nf'
include { CLEAN_OUTPUT_PATHS } from './modules/clean_output_paths.nf'
include { PACKAGE_PROCESSING_INFO } from './modules/package_processing_info.nf'
include { GENERATE_RAW_PROTOCOL } from './modules/generate_protocols.nf'
include { GENERATE_README } from './modules/generate_readme.nf'
include { GENERATE_MD5SUMS } from './modules/generate_md5sums.nf'
include { VALIDATE_PROCESSING } from './modules/validate_processing.nf'
include { UPDATE_ASSAY_TABLE } from './modules/update_assay_table.nf'


workflow {
    main:
    def processed_dir = file(params.outdir)
    if( !processed_dir.exists() ) {
            error "Output directory '${processed_dir}' does not exist. Make sure the main workflow is run first."
    }

    ch_runsheet = channel.fromPath("${processed_dir}/Metadata/*_runsheet.csv")

    ch_multiqc_data = channel.fromPath("${processed_dir}/*/MultiQC_Reports/*_data", type: 'dir')

    ch_unpurged_files = channel.fromPath("${processed_dir}/Aligned_Sequence_Data/Alignment_Reports/*report.txt")

    ch_processing_info = channel.fromPath("$launchDir/processing_scripts", type: 'dir')

    ch_software_versions = channel.fromPath("${processed_dir}/GeneLab/software_versions_*.md")

    ch_isa_zip = channel.fromPath("${processed_dir}/Metadata/*-ISA.zip")
    

    PARSE_RUNSHEET( ch_runsheet )
    samples = PARSE_RUNSHEET.out.samples
    runsheet_path = PARSE_RUNSHEET.out.runsheet

    // Extract metadata from the first sample and set it as a channel
    samples | first
            | map { meta, reads -> meta }
            | set { ch_meta }

    // Extract assay suffix from metadata to be used in processes where meta is not passed
    ch_meta | map { meta -> meta.assay_suffix }
            | set { assay_suffix }

    ch_multiqc_data_mapped = ch_multiqc_data.map { file -> [file, file.parent] }
    CLEAN_MULTIQC_PATHS( ch_multiqc_data_mapped )

    ch_unpurged_files_mapped = ch_unpurged_files.map { file -> [file, file.parent] }
    CLEAN_OUTPUT_PATHS( ch_unpurged_files_mapped )
    
    PACKAGE_PROCESSING_INFO( ch_processing_info, assay_suffix )

    // Run only when params.post_processing.include_raw_data is true, otherwise skip
    GENERATE_RAW_PROTOCOL( ch_meta, ch_software_versions )
    
    GENERATE_README( ch_meta )

    ch_ready = GENERATE_README.out.readme                    // README must exist
                .combine( CLEAN_MULTIQC_PATHS.out.collect() )   // Md5sum for _data.zip files from cleaned multiqc paths 
                .combine( CLEAN_OUTPUT_PATHS.out.collect() )   // Md5sum for cleaned alignment report txt files
                .map { true }
    
    GENERATE_MD5SUMS( assay_suffix, channel.value(processed_dir), ch_ready )

    VALIDATE_PROCESSING( ch_meta, runsheet_path, channel.value(processed_dir), ch_ready )

    ch_raw_multiqc_zip = CLEAN_MULTIQC_PATHS.out
                            .flatten()
                            .filter { zip -> zip.name.contains('raw_multiqc') }

    UPDATE_ASSAY_TABLE( ch_meta, runsheet_path, ch_isa_zip, ch_raw_multiqc_zip )  

    publish:
    processing_info = PACKAGE_PROCESSING_INFO.out.zip
    raw_protocol = GENERATE_RAW_PROTOCOL.out.raw_protocol
    readme = GENERATE_README.out.readme
    raw_md5sum = GENERATE_MD5SUMS.out.raw_md5sum
    processed_md5sum = GENERATE_MD5SUMS.out.processed_md5sum
    validation_log = VALIDATE_PROCESSING.out.validation_log
    assay_table = UPDATE_ASSAY_TABLE.out.assay_table
}

output {
    processing_info { 
        path "GeneLab" }

    raw_protocol {
        path "GeneLab" }

    readme {
        path "GeneLab" }

    raw_md5sum {
        path "GeneLab"
    }

    processed_md5sum {
        path "GeneLab" }

    validation_log {
        path "GeneLab" }

    assay_table {
        path "GeneLab/updated_curation_tables"
    }

}