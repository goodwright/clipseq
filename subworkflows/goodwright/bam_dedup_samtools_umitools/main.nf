/*
* UMIcollapse, index BAM file and run samtools stats, flagstat and idxstats
*/

include { UMICOLLAPSE        } from '../../../modules/nf-core/umicollapse/main'   
include { SAMTOOLS_INDEX     } from '../../../modules/nf-core/samtools/index/main'

workflow BAM_DEDUP_SAMTOOLS_UMITOOLS {
    take:
    bam_bai // channel: [ val(meta), [ bam ], [ bai/csi ] ]

    main:
    ch_versions = Channel.empty()

    /*
    * MODULE: UMI-tools collapse
    */
    UMICOLLAPSE ( 
        bam_bai,
        'bam' 
    )
    ch_versions = ch_versions.mix(UMICOLLAPSE.out.versions)

    /*
    * MODULE: Index BAM file
    */
    SAMTOOLS_INDEX ( 
        UMICOLLAPSE.out.bam 
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    UMICOLLAPSE.out.bam
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


    emit:
    bam      = UMICOLLAPSE.out.bam       // channel: [ val(meta), [ bam ] ]
    umi_log  = UMICOLLAPSE.out.log       // channel: [ val(meta), [ log ] ]
    bai      = SAMTOOLS_INDEX.out.bai          // channel: [ val(meta), [ bai ] ]
    csi      = SAMTOOLS_INDEX.out.csi          // channel: [ val(meta), [ csi ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}
