include { PARSE_RUNSHEET } from './parse_runsheet.nf'
include { FETCH_ISA } from '../modules/fetch_isa.nf'
include { ISA_TO_RUNSHEET } from '../modules/isa_to_runsheet.nf'
include { GET_ACCESSIONS } from '../modules/get_accessions.nf'
include { STAGE_RAW_READS } from './stage_raw_reads.nf'
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'
/**
 * STAGE_ANALYSIS
 * 
 * This subworkflow handles the initial setup of the Methylseq analysis:
 * 1. Fetches accessions if needed
 * 2. Obtains or creates the runsheet
 * 3. Parses the runsheet and stages raw reads
 */
workflow STAGE_ANALYSIS {
    take:
        dp_tools_plugin
        accession
        isa_archive_path
        runsheet_path
        api_url

    main:
        // Parse accession
        channel.empty() | set { osd_accession }
        channel.empty() | set { glds_accession }
        
        if ( accession ) {
            GET_ACCESSIONS( accession, api_url )
            osd_accession = GET_ACCESSIONS.out.accessions_txt.map { it.readLines()[0].trim() }
            glds_accession = GET_ACCESSIONS.out.accessions_txt.map { it.readLines()[1].trim() }
        }

        channel.empty() | set { isa_archive }
        channel.empty() | set { dp_tools_version }
        if ( runsheet_path == null ) { // if runsheet_path is not provided, set it up from ISA input
            if ( isa_archive_path == null ) { // if isa_archive_path is not provided, fetch the ISA
                FETCH_ISA( osd_accession, glds_accession )
                isa_archive = FETCH_ISA.out.isa_archive
            } else {
                // isa_archive_path is already a channel, use it directly
                isa_archive = isa_archive_path
            }
            ISA_TO_RUNSHEET( osd_accession, glds_accession, isa_archive, dp_tools_plugin )
            runsheet_path = ISA_TO_RUNSHEET.out.runsheet
            dp_tools_version = ISA_TO_RUNSHEET.out.version
        }

        // Validate input parameters and runsheet
        validateParameters()

        PARSE_RUNSHEET( runsheet_path )
        samples = PARSE_RUNSHEET.out.samples
        runsheet_path = PARSE_RUNSHEET.out.runsheet

        // Stage the full or truncated raw reads
        STAGE_RAW_READS( samples )
        raw_reads = STAGE_RAW_READS.out.raw_reads
        
    emit:
        samples         = samples
        raw_reads       = raw_reads
        runsheet_path   = runsheet_path
        isa_archive     = isa_archive
        osd_accession   = osd_accession
        glds_accession  = glds_accession
        dp_tools_version = dp_tools_version
} 