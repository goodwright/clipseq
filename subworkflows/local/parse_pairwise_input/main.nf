/*
 * Get sample id from pairwise samplesheet and check if compatible with fastq samplesheet
 * Then create channel of pairwise ids for bam file channel reorganisation before MACS3 peak calling
 * Subworkflow goal is to create pairwise channels of sample id that match metadata
 */

include { PAIRWISE_SAMPLESHEET_CHECK } from '../../../modules/local/whole_read_analysis/pairwise_samplesheet_check/main.nf'

workflow PARSE_PAIRWISE_INPUT {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    ch_versions = Channel.empty()

    /*
    * MODULE: Check the samplesheet for errors
    */
    PAIRWISE_SAMPLESHEET_CHECK (
        samplesheet
    )

    /*
    * CHANNEL: Split the validated pairwise samplesheet into a channel
    */
    ch_pairwise_samples = PAIRWISE_SAMPLESHEET_CHECK.out.csv
        .splitCsv ( header:true, sep:"," )
        .map { row -> [row.id, row.control_id] }

    emit:
    pairwise = ch_pairwise_samples // channel: [ id, control_id ]
    versions = ch_versions  // channel: [ versions.yml ]
}
