/* Split bam files into plus and minus strand bam files for MACS3 peak calling
*  Also convert the bam file to bed format and calculate (normalised) coverage tracks in bedgraph format

* Update: The new workflow generates bigwig files from positive and negative strand bedgraph files
*/

/*
* MODULES
*/
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_PLUS          } from '../../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_MINUS         } from '../../../modules/nf-core/samtools/view/main'

include { BEDTOOLS_BAMTOBED                            } from '../../../modules/nf-core/bedtools/bamtobed/main.nf'
include { BEDTOOLS_GENOMECOV as BEDTOOLS_GENOMECOV_POS } from '../../../modules/nf-core/bedtools/genomecov/main.nf'
include { BEDTOOLS_GENOMECOV as BEDTOOLS_GENOMECOV_NEG } from '../../../modules/nf-core/bedtools/genomecov/main.nf'

include { UCSC_BEDGRAPHTOBIGWIG as BIGWIG_POS          } from '../../../modules/nf-core/ucsc/bedgraphtobigwig/main.nf'
include { UCSC_BEDGRAPHTOBIGWIG as BIGWIG_NEG          } from '../../../modules/nf-core/ucsc/bedgraphtobigwig/main.nf'

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
    * MODULE: Report depth on the pos/neg strand, using whole interval as bedgraph output
    */
    BEDTOOLS_GENOMECOV_POS (
        BEDTOOLS_BAMTOBED.out.bed.map{ [ it[0], it[1], 1 ] },
        fai,
        'pos.bedgraph'
    )
    ch_versions = ch_versions.mix(BEDTOOLS_GENOMECOV_POS.out.versions)

    BEDTOOLS_GENOMECOV_NEG (
        BEDTOOLS_BAMTOBED.out.bed.map{ [ it[0], it[1], 1 ] },
        fai,
        'neg.bedgraph'
    )

    /*
    * MODULE: Convert BED file to bigwig format
    */
    BIGWIG_POS (
        BEDTOOLS_GENOMECOV_POS.out.genomecov,
        params.chrom_sizes
    )
    ch_versions = ch_versions.mix(BIGWIG_POS.out.versions)

    BIGWIG_NEG (
        BEDTOOLS_GENOMECOV_NEG.out.genomecov,
        params.chrom_sizes
    )

    emit:
    bam_plus     = SAMTOOLS_VIEW_PLUS.out.bam       // channel: [ val(meta), [ bam ] ]
    bam_minus    = SAMTOOLS_VIEW_MINUS.out.bam      // channel: [ val(meta), [ bam ] ]

    bed           = BEDTOOLS_BAMTOBED.out.bed         // channel: [ val(meta), [ bed ] ]
    versions      = ch_versions                     // channel: [ versions.yml ]
}
