// Include the COPY_READS module
include { COPY_READS } from '../modules/copy_reads.nf'

workflow STAGE_RAW_READS {
    take:
        ch_samples

    main:

        if (params.truncate_to) {

            // Download truncated raw reads
            ch_samples | map { it -> it[0].paired_end ? [it[0], it[1][0], it[1][1]] : [it[0], it[1][0]]}
                | branch {
                paired: it.size() == 3
                single: it.size() == 2
                }
                | set{ ch_raw_read_pointers }
            
            
            // TO DO: Move the two splitFastq truncation steps into processes that run in parallel. Low priority since this is just for debugging.

            // PAIRED END
            // Only difference is the splitFastq arg 'pe'
            ch_raw_read_pointers.paired | splitFastq(pe: true, decompress: true, compress: true, limit: params.truncate_to, by: params.truncate_to, file: true)
                                        | map { it -> [ it[0], [ it[1], it[2] ] ]}
                                        | set { ch_raw_reads }
            // SINGLE END
            // Only difference is the splitFastq arg 'pe'
            ch_raw_read_pointers.single | splitFastq(decompress: true, compress: true, limit: params.truncate_to, by: params.truncate_to, file: true)
                                        | map { it -> [ it[0], [ it[1] ] ]}
                                        | mix( ch_raw_reads )
                                        | set { ch_raw_reads }

            // Moves the truncated files to expected raw read locations as per samplesheet
            ch_raw_reads | COPY_READS

        } else {
            ch_samples | map { it -> it[0].paired_end ? [it[0], [ it[1][0], it[1][1] ]] : [it[0], [it[1][0]]]}
                | set { ch_raw_reads }

            ch_raw_reads | COPY_READS
        } 

    emit:
        raw_reads = COPY_READS.out.raw_reads
}
