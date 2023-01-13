//
// Read QC, read trimming and post trim QC
//

include { FASTQC     } from '../../../modules/nf-core/fastqc/main'
include { TRIMGALORE } from '../../../modules/nf-core/trimgalore/main'

workflow FASTQC_TRIMGALORE {
    take:
    fastq         // channel: [ val(meta), [ fastq ] ]
    skip_fastqc   // boolean: true/false
    skip_trimming // boolean: true/false

    main:
    ch_versions = Channel.empty()

    /*
    * MODULE: Run pre-trim fastqc
    */
    ch_fastqc_html = Channel.empty()
    ch_fastqc_zip  = Channel.empty()
    if (!skip_fastqc) {
        FASTQC (
            fastq
        )
        ch_fastqc_zip  = FASTQC.out.zip
        ch_fastqc_html = FASTQC.out.html
        ch_versions    = ch_versions.mix(FASTQC.out.versions)
    }

    /*
    * MODULE: Trim reads and run post-trim fastqc
    */
    ch_fastqc_trim_html = Channel.empty()
    ch_fastqc_trim_zip  = Channel.empty()
    ch_trim_log         = Channel.empty()
    ch_output_reads     = fastq
    if (!skip_trimming) {
        TRIMGALORE (
            fastq
        )
        ch_output_reads     = TRIMGALORE.out.reads
        ch_fastqc_trim_html = TRIMGALORE.out.html
        ch_fastqc_trim_zip  = TRIMGALORE.out.zip
        ch_trim_log         = TRIMGALORE.out.log
        ch_versions         = ch_versions.mix(TRIMGALORE.out.versions)
    }

    emit:
    fastq            = ch_output_reads     // channel: [ val(meta), [ reads ] ]
    versions         = ch_versions         // channel: [ versions.yml ]
    fastqc_html      = ch_fastqc_html      // channel: [ val(meta), [ html ] ]
    fastqc_zip       = ch_fastqc_zip       // channel: [ val(meta), [ zip ] ]
    fastqc_trim_html = ch_fastqc_trim_html // channel: [ val(meta), [ html ] ]
    fastqc_trim_zip  = ch_fastqc_trim_zip  // channel: [ val(meta), [ zip ] ]
    trim_log         = ch_trim_log         // channel: [ val(meta), [ txt ] ]
}
