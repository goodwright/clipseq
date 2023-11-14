//
// Filters fastq files for unwanted alignments before aligning to genome and transcriptome using STAR
// BAM alignments are then sorted, indexed and stats are calculated
//

/*
* MODULES
*/ 
include { BOWTIE_ALIGN                                } from '../../../modules/nf-core/bowtie/align/main.nf'
include { STAR_ALIGN                                  } from '../../../modules/nf-core/star/align/main.nf'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_TRANSCRIPT   } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_SMRNA        } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_SMRNA      } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_TRANSCRIPT } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_GENOME     } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_GENOME       } from '../../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_TRANSCRIPT   } from '../../../modules/nf-core/samtools/view/main'  

workflow RNA_ALIGN {
    take:
    fastq      // channel: [ val(meta), [ fastq ] ]
    bt_index   // channel: [ val(meta), index ]
    star_index // channel: [ val(meta), index ]
    gtf        // channel: [ val(meta), gtf ]
    fasta      // channel: [ val(meta), fasta/fa ]

    main:
    ch_versions = Channel.empty()

    /*
    * MODULE: Align reads to smrna genome
    */
    BOWTIE_ALIGN (
        fastq,
        bt_index.collect{it[1]}
    )
    ch_versions = ch_versions.mix(BOWTIE_ALIGN.out.versions)

    SAMTOOLS_SORT_SMRNA ( BOWTIE_ALIGN.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT_SMRNA.out.versions)

    SAMTOOLS_INDEX_SMRNA ( SAMTOOLS_SORT_SMRNA.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX_SMRNA.out.versions)

    /*
    * MODULE: Align reads that did not align to the smrna genome to the primary genome
    */
    STAR_ALIGN (
        BOWTIE_ALIGN.out.fastq,
        star_index.collect{ [it[0], it[1]]},
        gtf.collect{ [it[0], it[1]]},
        false,
        '',
        ''
    )
    ch_versions = ch_versions.mix(STAR_ALIGN.out.versions)

    /*
    * MODULE: Index genome-level BAM file 
    */
    SAMTOOLS_INDEX_GENOME ( STAR_ALIGN.out.bam_sorted )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX_GENOME.out.versions.first())

    /*
    * MODULE: Index transcript-level BAM file 
    */
    SAMTOOLS_SORT_TRANSCRIPT ( STAR_ALIGN.out.bam_transcript )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT_TRANSCRIPT.out.versions.first())

    SAMTOOLS_INDEX_TRANSCRIPT ( SAMTOOLS_SORT_TRANSCRIPT.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX_TRANSCRIPT.out.versions.first())


    /*
    * CHANNEL: Join bam and bai files
    */
    ch_bam_bai = STAR_ALIGN.out.bam_sorted
        .join(SAMTOOLS_INDEX_GENOME.out.bai, by: [0], remainder: true)
        .join(SAMTOOLS_INDEX_GENOME.out.csi, by: [0], remainder: true)
        .map {
            meta, bam, bai, csi ->
                if (bai) {
                    [ meta, bam, bai ]
                } else {
                    [ meta, bam, csi ]
                }
        }

    /*
    * CHANNEL: Join bam and bai files
    */
    ch_transcript_bam_bai = SAMTOOLS_SORT_TRANSCRIPT.out.bam
        .join(SAMTOOLS_INDEX_TRANSCRIPT.out.bai, by: [0], remainder: true)
        .join(SAMTOOLS_INDEX_TRANSCRIPT.out.csi, by: [0], remainder: true)
        .map {
            meta, bam, bai, csi ->
                if (bai) {
                    [ meta, bam, bai ]
                } else {
                    [ meta, bam, csi ]
                }
        }


    /*
    * CHANNEL: Filter for uniquely mapping reads for downstream analysis; samtools view -b -q 5 -o output.bam alignments.bam
    */
    SAMTOOLS_VIEW_GENOME (
        ch_bam_bai,
        [[],[]],
        []
    )

    SAMTOOLS_VIEW_TRANSCRIPT (
        ch_transcript_bam_bai,
        [[],[]],
        []
    )

    emit:
    bt_bam              = BOWTIE_ALIGN.out.bam                            // channel: [ val(meta), [ bam ] ]
    bt_log              = BOWTIE_ALIGN.out.log                            // channel: [ val(meta), [ txt ] ]
    star_log            = STAR_ALIGN.out.log                              // channel: [ val(meta), [ txt ] ]
    star_log_final      = STAR_ALIGN.out.log_final                        // channel: [ val(meta), [ txt ] ]
    genome_unique_bam   = SAMTOOLS_VIEW_GENOME.out.bam                    // channel: [ val(meta), [ bam ] ]
    genome_unique_bai   = SAMTOOLS_VIEW_GENOME.out.bai                    // channel: [ val(meta), [ bai ] ]
    genome_multi_bam    = STAR_ALIGN.out.bam_sorted                       // channel: [ val(meta), [ bam ] ]
    genome_multi_bai    = SAMTOOLS_INDEX_GENOME.out.bai                   // channel: [ val(meta), [ bai ] ]
    transcript_bam      = SAMTOOLS_VIEW_TRANSCRIPT.out.bam                // channel: [ val(meta), [ bam ] ]
    transcript_bai      = SAMTOOLS_VIEW_TRANSCRIPT.out.bai                // channel: [ val(meta), [ bai ] ]
    smrna_bam           = SAMTOOLS_SORT_SMRNA.out.bam                     // channel: [ val(meta), [ bam ] ]
    smrna_bai           = SAMTOOLS_INDEX_SMRNA.out.bai                    // channel: [ val(meta), [ bai ] ]
    versions            = ch_versions                                     // channel: [ versions.yml ]
}
