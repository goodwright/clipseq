//
// Filters fastq files for unwanted alignments before aligning to genome and transcriptome using STAR
// BAM alignments are then sorted, indexed and stats are calculated
//

/*
* MODULES
*/
include { BOWTIE_ALIGN                                } from '../../../modules/nf-core/bowtie/align/main.nf'
include { STAR_ALIGN                                  } from '../../../modules/nf-core/star/align/main.nf'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_GENOME     } from '../../../modules/nf-core/samtools/index/main.nf'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_TRANSCRIPT } from '../../../modules/nf-core/samtools/index/main.nf'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_TRANSCRIPT   } from '../../../modules/nf-core/samtools/sort/main.nf'


workflow RNA_ALIGN {
    take:
    fastq      // channel: [ val(meta), [ fastq ] ]
    bt2_index  // channel: [ index ]
    star_index // channel: [ index ]
    gtf        // channel: [ gtf ]
    fasta      // channel: [ fasta/fa ]

    main:
    ch_versions = Channel.empty()

    /*
    * MODULE: Align reads to smrna genome
    */
    BOWTIE_ALIGN (
        fastq,
        bt2_index
    )
    ch_versions = ch_versions.mix(BOWTIE_ALIGN.out.versions)

    /*
    * MODULE: Align reads that did not align to the smrna genome to the primary genome
    */
    STAR_ALIGN (
        BOWTIE_ALIGN.out.fastq,
        star_index,
        gtf,
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

    SAMTOOLS_SORT_TRANSCRIPT ( STAR_ALIGN.out.bam_transcript )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT_TRANSCRIPT.out.versions.first())

    SAMTOOLS_INDEX_TRANSCRIPT ( SAMTOOLS_SORT_TRANSCRIPT.out.bam )

    emit:
    bt_bam              = BOWTIE_ALIGN.out.bam                            // channel: [ val(meta), [ bam ] ]
    bt_log              = BOWTIE_ALIGN.out.log                            // channel: [ val(meta), [ txt ] ]
    star_bam            = STAR_ALIGN.out.bam_sorted                       // channel: [ val(meta), [ bam ] ]
    star_bam_transcript = STAR_ALIGN.out.bam_transcript                   // channel: [ val(meta), [ bam ] ]
    star_log            = STAR_ALIGN.out.log                              // channel: [ val(meta), [ txt ] ]
    star_log_final      = STAR_ALIGN.out.log_final                        // channel: [ val(meta), [ txt ] ]
    genome_bam          = STAR_ALIGN.out.bam_sorted                       // channel: [ val(meta), [ bam ] ]
    genome_bai          = SAMTOOLS_INDEX_GENOME.out.bai                          // channel: [ val(meta), [ bai ] ]
    transcript_bam      = SAMTOOLS_SORT_TRANSCRIPT.out.bam      // channel: [ val(meta), [ bam ] ]
    transcript_bai      = SAMTOOLS_INDEX_TRANSCRIPT.out.bai      // channel: [ val(meta), [ bai ] ]
    versions            = ch_versions                                     // channel: [ versions.yml ]
}
