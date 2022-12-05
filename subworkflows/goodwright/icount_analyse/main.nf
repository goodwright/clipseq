//
// Runs suite of icount tools on an input crosslinks BED file
//

/*
* MODULES
*/
include { ICOUNT_SUMMARY          } from '../../../modules/goodwright/icount/summary/main.nf'
include { ICOUNT_RNAMAPS          } from '../../../modules/goodwright/icount/rnamaps/main.nf'
include { ICOUNT_SIGXLS           } from '../../../modules/goodwright/icount/sigxls/main.nf'
include { ICOUNT_PEAKS            } from '../../../modules/goodwright/icount/peaks/main.nf'
include { GUNZIP as GUNZIP_SIGXLS } from '../../../modules/nf-core/gunzip/main.nf'
include { GUNZIP as GUNZIP_PEAKS  } from '../../../modules/nf-core/gunzip/main.nf'

workflow ICOUNT_ANALYSE {
    take:
    bed             // channel: [ val(meta), [ bed ] ]
    gtf_seg         // channel: [ [ gtf ] ]
    gtf_resolved    // channel: [ [ gtf.gz ] ]
    run_peakcalling // val: boolean

    main:
    ch_versions = Channel.empty()

    /*
    * MODULE: Run iCount summary
    */
    ICOUNT_SUMMARY (
        bed,
        gtf_resolved
    )
    ch_versions = ch_versions.mix(ICOUNT_SUMMARY.out.versions)

    /*
    * MODULE: Run iCount rnamaps
    */
    ICOUNT_RNAMAPS (
        bed,
        gtf_resolved
    )
    ch_versions = ch_versions.mix(ICOUNT_RNAMAPS.out.versions)

    ch_bed_sigxls = Channel.empty()
    ch_tsv_scores = Channel.empty()
    ch_bed_peaks  = Channel.empty()
    if (run_peakcalling) {
        /*
        * MODULE: Run iCount sigxls
        */
        ICOUNT_SIGXLS (
            bed,
            gtf_seg
        )
        ch_versions   = ch_versions.mix(ICOUNT_SIGXLS.out.versions)
        ch_tsv_scores = ICOUNT_SIGXLS.out.scores

        /*
        * CHANNEL: Create combined channel of input crosslinks and sigxls
        */
        ch_peaks_input = bed
            .map{ [ it[0].id, it[0], it[1] ] }
            .join( ICOUNT_SIGXLS.out.sigxls.map{ [ it[0].id, it[0], it[1] ] } )
            .map { [ it[1], it[2], it[4] ] }
        //EXAMPLE CHANNEL STRUCT: [ [id:test], BED(crosslinks), BED.GZ(sigxls) ]
        //ch_peaks_input | view 

        /*
        * MODULE: Run iCount peaks
        */
        ICOUNT_PEAKS (
            ch_peaks_input
        )
        ch_versions = ch_versions.mix(ICOUNT_PEAKS.out.versions)

        /*
        * MODULE: Decompress sigxls and peaks
        */
        GUNZIP_SIGXLS (
            ICOUNT_SIGXLS.out.sigxls
        )
        GUNZIP_PEAKS (
            ICOUNT_PEAKS.out.peaks
        )
        ch_versions   = ch_versions.mix(GUNZIP_SIGXLS.out.versions)
        ch_bed_sigxls = GUNZIP_SIGXLS.out.gunzip
        ch_bed_peaks  = GUNZIP_PEAKS.out.gunzip
    }

    emit:
    tsv_summary_type    = ICOUNT_SUMMARY.out.summary_type    // channel: [ val(meta), [ tsv ] ]
    tsv_summary_subtype = ICOUNT_SUMMARY.out.summary_subtype // channel: [ val(meta), [ tsv ] ]
    tsv_summary_gene    = ICOUNT_SUMMARY.out.summary_gene    // channel: [ val(meta), [ tsv ] ]
    tsv_rnamaps         = ICOUNT_RNAMAPS.out.tsv             // channel: [ val(meta), [ tsv ] ]
    bed_sigxls          = ch_bed_sigxls                      // channel: [ val(meta), [ bed ] ]
    tsv_scores          = ch_tsv_scores                      // channel: [ val(meta), [ tsv ] ]
    bed_peaks           = ch_bed_peaks                       // channel: [ val(meta), [ bed ] ]
    versions            = ch_versions                        // channel: [ versions.yml ]
}
