/* Split bam files into plus and minus strand bam files for MACS3 peak calling
*  Also convert the bam file to bed format and calculate (normalised) coverage tracks in bedgraph format

* Update: The new workflow generates genomecov bedgraph files to save disk space and time
*         For this we use genomecov -bg option to generate bedgraph files directly
*         The bed files are per-read bed files converted from bam files
*         The bedgraph files are genomecov files converted from bed files
*         The coverage tracks are bedgraph files scaled by total coverage using the same method as before:
*         - total = total + (length * abs_val_of_bedgraph_value_on_the_region)
*         - data value = (data value / total) * 1e6
*/

/*
* MODULES
*/
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_PLUS          } from '../../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_MINUS         } from '../../../modules/nf-core/samtools/view/main'

include { BEDTOOLS_BAMTOBED                            } from '../../../modules/nf-core/bedtools/bamtobed/main.nf'
include { BEDTOOLS_GENOMECOV as BEDTOOLS_GENOMECOV_POS } from '../../../modules/nf-core/bedtools/genomecov/main.nf'
include { BEDTOOLS_GENOMECOV as BEDTOOLS_GENOMECOV_NEG } from '../../../modules/nf-core/bedtools/genomecov/main.nf'
include { LINUX_COMMAND as SELECT_BED_POS              } from '../../../modules/goodwright/linux/command/main.nf'
include { LINUX_COMMAND as SELECT_BED_NEG              } from '../../../modules/goodwright/linux/command/main.nf'
include { LINUX_COMMAND as MERGE_AND_SORT              } from '../../../modules/goodwright/linux/command/main.nf'
include { LINUX_COMMAND as CROSSLINK_COVERAGE          } from '../../../modules/goodwright/linux/command/main.nf'
include { LINUX_COMMAND as CROSSLINK_NORMCOVERAGE      } from '../../../modules/goodwright/linux/command/main.nf'

workflow PREP_WHOLE_READ {
    take:
    bam // channel: [ val(meta), [ bam ] ]
    bai // channel: [ val(meta), [ bai ] ]
    fai // channel: [ fai ]

    main:
    ch_versions = Channel.empty()

    /*
    * MODULE: split BAM file into plus and minus strand files
    */
    ch_bam = bam.join(bai, by: [0], remainder: true)
        .map { meta, bam, bai -> [ meta, bam, bai ] }

    SAMTOOLS_VIEW_PLUS (
        ch_bam,
        [[],[]],
        []
    )
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW_PLUS.out.versions)

    SAMTOOLS_VIEW_MINUS (
        ch_bam,
        [[],[]],
        []
    )

    /*
    * MODULE: Convert input BAM file to BED format
    */
    BEDTOOLS_BAMTOBED (
        bam
    )
    ch_versions = ch_versions.mix(BEDTOOLS_BAMTOBED.out.versions)

    /*
    * MODULE: Report depth at each position on the pos strand, using whole interval
    */
    BEDTOOLS_GENOMECOV_POS (
        BEDTOOLS_BAMTOBED.out.bed.map{ [ it[0], it[1], 1 ] },
        fai,
        'pos.bedgraph'
    )
    ch_versions = ch_versions.mix(BEDTOOLS_GENOMECOV_POS.out.versions)

    /*
    * MODULE: Report depth at each position on the neg strand, using whole interval
    */
    BEDTOOLS_GENOMECOV_NEG (
        BEDTOOLS_BAMTOBED.out.bed.map{ [ it[0], it[1], 1 ] },
        fai,
        'neg.bedgraph'
    )

    /*
    * CHANNEL: Join POS/NEG files into one channel so they can be merged in the next module
    */
    ch_merge_and_sort_input = BEDTOOLS_GENOMECOV_POS.out.genomecov
        .map{ [ it[0].id, it[0], it[1] ] }
        .join( BEDTOOLS_GENOMECOV_NEG.out.genomecov.map{ [ it[0].id, it[0], it[1] ] } )
        .map { [ it[1], [ it[2], it[4] ] ] }
    //EXAMPLE CHANNEL STRUCT: [ [id:test], [ BED(pos), BED(neg) ] ]
    //ch_merge_and_sort_input | view 

    /*
    * MODULE: Select columns in BED file using AWK
    */
    MERGE_AND_SORT (
        ch_merge_and_sort_input,
        [],
        false
    )

    /*
    * MODULE: Create normalised coverage track using AWK
    */
    CROSSLINK_NORMCOVERAGE (
        MERGE_AND_SORT.out.file,
        [],
        true
    )

    emit:
    bam_plus     = SAMTOOLS_VIEW_PLUS.out.bam       // channel: [ val(meta), [ bam ] ]
    bam_minus    = SAMTOOLS_VIEW_MINUS.out.bam      // channel: [ val(meta), [ bam ] ]

    bed           = BEDTOOLS_BAMTOBED.out.bed         // channel: [ val(meta), [ bed ] ]
    coverage      = MERGE_AND_SORT.out.file     // channel: [ val(meta), [ bedgraph ] ]
    coverage_norm = CROSSLINK_NORMCOVERAGE.out.file // channel: [ val(meta), [ bedgraph ] ]
    versions      = ch_versions                     // channel: [ versions.yml ]
}
