/*
* UMI-tools collapse, index BAM file and run samtools stats, flagstat and idxstats
*/

include { UMITOOLS_COLLAPSE  } from '../../../modules/goodwright/umitools/collapse/main'
include { SAMTOOLS_INDEX     } from '../../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS } from '../../nf-core/bam_stats_samtools/main'

workflow BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS {
    take:
    bam_bai // channel: [ val(meta), [ bam ], [ bai/csi ] ]

    main:
    ch_versions = Channel.empty()

    /*
    * MODULE: UMI-tools collapse
    */
    UMITOOLS_COLLAPSE ( 
        bam_bai 
    )
    ch_versions = ch_versions.mix(UMITOOLS_COLLAPSE.out.versions)

    /*
    * MODULE: Index BAM file
    */
    SAMTOOLS_INDEX ( 
        UMITOOLS_COLLAPSE.out.bam 
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    UMITOOLS_COLLAPSE.out.bam
        .join(SAMTOOLS_INDEX.out.bai, by: [0], remainder: true)
        .join(SAMTOOLS_INDEX.out.csi, by: [0], remainder: true)
        .map {
            meta, bam, bai, csi ->
                if (bai) {
                    [ meta, bam, bai ]
                } else {
                    [ meta, bam, csi ]
                }
        }
        .set { ch_bam_bai }

    /*
    * SUBWORKFLOW: Calc stats on new bam
    */
    BAM_STATS_SAMTOOLS ( 
        ch_bam_bai, [] 
    )
    ch_versions = ch_versions.mix(BAM_STATS_SAMTOOLS.out.versions)

    emit:
    bam      = UMITOOLS_COLLAPSE.out.bam       // channel: [ val(meta), [ bam ] ]
    umi_log  = UMITOOLS_COLLAPSE.out.log       // channel: [ val(meta), [ log ] ]
    bai      = SAMTOOLS_INDEX.out.bai          // channel: [ val(meta), [ bai ] ]
    csi      = SAMTOOLS_INDEX.out.csi          // channel: [ val(meta), [ csi ] ]
    stats    = BAM_STATS_SAMTOOLS.out.stats    // channel: [ val(meta), [ stats ] ]
    flagstat = BAM_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats = BAM_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}
