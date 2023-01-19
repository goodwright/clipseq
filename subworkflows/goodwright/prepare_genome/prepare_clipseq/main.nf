//
// Prepare all genome files for running the clipseq analysis pipeline
//

/*
* MODULES
*/
include { CLIPSEQ_FIND_LONGEST_TRANSCRIPT                                        } from '../../../../modules/goodwright/clipseq/find_longest_transcript/main.nf'
include { CLIPSEQ_FILTER_GTF                                                     } from '../../../../modules/goodwright/clipseq/filter_gtf/main.nf'
include { ICOUNT_SEGMENT as ICOUNT_SEG_GTF                                       } from '../../../../modules/goodwright/icount/segment/main.nf'
include { ICOUNT_SEGMENT as ICOUNT_SEG_FILTGTF                                   } from '../../../../modules/goodwright/icount/segment/main.nf'
include { CLIPSEQ_RESOLVE_UNANNOTATED as RESOLVE_UNANNOTATED                     } from '../../../../modules/goodwright/clipseq/resolve_unannotated/main.nf'
include { CLIPSEQ_RESOLVE_UNANNOTATED as RESOLVE_UNANNOTATED_GENIC_OTHER         } from '../../../../modules/goodwright/clipseq/resolve_unannotated/main.nf'
include { CLIPSEQ_RESOLVE_UNANNOTATED as RESOLVE_UNANNOTATED_REGIONS             } from '../../../../modules/goodwright/clipseq/resolve_unannotated/main.nf'
include { CLIPSEQ_RESOLVE_UNANNOTATED as RESOLVE_UNANNOTATED_GENIC_OTHER_REGIONS } from '../../../../modules/goodwright/clipseq/resolve_unannotated/main.nf'

/*
* SUBWORKFLOWS
*/
include { PREPARE_REF as PREPARE_PRIMARY_GENOME    } from '../prepare_ref/main.nf'
include { PREPARE_REF as PREPARE_SMRNA_GENOME      } from '../prepare_ref/main.nf'
include { PREPARE_ALINGER as PREPARE_PRIMARY_INDEX } from '../prepare_aligner/main.nf'
include { PREPARE_ALINGER as PREPARE_SMRNA_INDEX   } from '../prepare_aligner/main.nf'

workflow PREPARE_CLIPSEQ {
    take:
    fasta                            // channel: [ fasta ]
    smrna_fasta                      // channel: [ fasta ]
    gtf                              // channel: [ gtf ]
    genome_index_path                // channel: [ folder/tar.gz ]
    smrna_index_path                 // channel: [ folder/tar.gz ]
    fasta_fai                        // channel: [ fasta.fai ]                                            
    filtered_gtf                     // channel: [ gtf ]                      
    chrom_sizes                      // channel: [ tsv ]                     
    smrna_fasta_fai                  // channel: [ fasta.fai ]                         
    smrna_chrom_sizes                // channel: [ tsv ]                           
    longest_transcript               // channel: [ txt ]                            
    seg_gtf                          // channel: [ gtf ]                 
    seg_filt_gtf                     // channel: [ gtf ]                      
    seg_resolved_gtf                 // channel: [ gtf ]                          
    seg_resolved_gtf_genic           // channel: [ gtf ]                                
    regions_gtf                      // channel: [ gtf ]                     
    regions_filt_gtf                 // channel: [ gtf ]                          
    regions_resolved_gtf             // channel: [ gtf ]                              
    regions_resolved_gtf_genic       // channel: [ gtf ]                                    
    longest_transcript_fai           // channel: [ fasta.fai ]                                

    main:
    ch_versions = Channel.empty()

    /*
    * SUBWORKFLOW: Uncompress and prepare main genome files
    */
    if (params.fasta_fai && params.chrom_sizes){
        ch_fasta       = params.fasta
        ch_fasta_fai   = params.fasta_fai
        ch_gtf         = params.gtf
        ch_chrom_sizes = params.chrom_sizes
    } else {
        PREPARE_PRIMARY_GENOME (
            fasta,
            gtf,
            [],
            []
        )
        ch_fasta       = PREPARE_PRIMARY_GENOME.out.fasta
        ch_fasta_fai   = PREPARE_PRIMARY_GENOME.out.fasta_fai
        ch_gtf         = PREPARE_PRIMARY_GENOME.out.gtf
        ch_chrom_sizes = PREPARE_PRIMARY_GENOME.out.chrom_sizes
        ch_versions    = ch_versions.mix(PREPARE_PRIMARY_GENOME.out.versions)
    }

    /*
    * SUBWORKFLOW: Uncompress and prepare smrna genome files
    */
    if (params.smrna_fasta_fai && params.smrna_chrom_sizes){
        ch_smrna_fasta       = params.smrna_fasta
        ch_smrna_fasta_fai   = params.smrna_fasta_fai
        ch_smrna_chrom_sizes = params.smrna_chrom_sizes
    } else {
        PREPARE_SMRNA_GENOME (
            smrna_fasta,
            [],
            [],
            []
        )
        ch_smrna_fasta       = PREPARE_SMRNA_GENOME.out.fasta
        ch_smrna_fasta_fai   = PREPARE_SMRNA_GENOME.out.fasta_fai
        ch_smrna_chrom_sizes = PREPARE_SMRNA_GENOME.out.chrom_sizes
        ch_versions          = ch_versions.mix(PREPARE_SMRNA_GENOME.out.versions)
    }

    /*
    * MODULE: Find the longest transcript from the primary genome
    */
    if (params.longest_transcript && params.longest_transcript_fai){
        ch_longest_transcript     = params.longest_transcript
        ch_longest_transcript_fai = params.longest_transcript_fai
    } else {
        CLIPSEQ_FIND_LONGEST_TRANSCRIPT (
            ch_gtf
        )
        ch_longest_transcript     = CLIPSEQ_FIND_LONGEST_TRANSCRIPT.out.longest_transcript
        ch_longest_transcript_fai = CLIPSEQ_FIND_LONGEST_TRANSCRIPT.out.longest_transcript_fai
        ch_versions               = ch_versions.mix(CLIPSEQ_FIND_LONGEST_TRANSCRIPT.out.versions)
    }

    /*
    * MODULE: Filter the GTF file
    */
    if (params.filtered_gtf){
    ch_filt_gtf = params.filtered_gtf
    } else {
    CLIPSEQ_FILTER_GTF (
        ch_gtf
    )
    ch_filt_gtf = CLIPSEQ_FILTER_GTF.out.gtf
    ch_versions = ch_versions.mix(CLIPSEQ_FILTER_GTF.out.versions)
    }

    /*
    * SUBWORKFLOW: Prepare STAR index for primary genome
    */
    PREPARE_PRIMARY_INDEX (
        ["star"],
        ch_fasta,
        ch_gtf,
        [],
        genome_index_path
    )
    ch_genome_index = PREPARE_PRIMARY_INDEX.out.star_index
    ch_versions     = ch_versions.mix(PREPARE_PRIMARY_INDEX.out.versions)

    /*
    * SUBWORKFLOW: Prepare BT2 index for smrna genome
    */
    PREPARE_SMRNA_INDEX (
        ["bowtie"],
        ch_smrna_fasta,
        [],
        smrna_index_path,
        []
    )
    ch_smrna_index = PREPARE_SMRNA_INDEX.out.bt2_index
    ch_versions    = ch_versions.mix(PREPARE_SMRNA_INDEX.out.versions)

    /*
    * MODULE: Segment GTF file using icount
    */
    if (params.seg_gtf && params.regions_gtf){
    ch_seg_gtf     = params.seg_gtf
	ch_regions_gtf = params.regions_gtf
    } else {
    ICOUNT_SEG_GTF (
        ch_gtf,
        ch_fasta_fai.map{ it[1] }
    )
    ch_seg_gtf     = ICOUNT_SEG_GTF.out.gtf
	ch_regions_gtf = ICOUNT_SEG_GTF.out.regions
    ch_versions    = ch_versions.mix(ICOUNT_SEG_GTF.out.versions)
    }

    /*
    * MODULE: Segment the filtered GTF file using icount
    */
    if (params.seg_filt_gtf && params.regions_filt_gtf){
    ch_seg_filt_gtf     = params.seg_filt_gtf
    ch_regions_filt_gtf = params.regions_filt_gtf
    } else {
    ICOUNT_SEG_FILTGTF (
        ch_filt_gtf,
        ch_fasta_fai.map{ it[1] }
    )
    ch_seg_filt_gtf     = ICOUNT_SEG_FILTGTF.out.gtf
    ch_regions_filt_gtf = ICOUNT_SEG_FILTGTF.out.regions
    }

    /*
    * MODULE: Resolve the GTF regions that iCount did not annotate
    */
    if (params.seg_resolved_gtf){
    ch_seg_resolved_gtf = params.seg_resolved_gtf
    } else {
    RESOLVE_UNANNOTATED (
        ICOUNT_SEG_GTF.out.gtf.map{ it[1] },
        ICOUNT_SEG_FILTGTF.out.gtf.map{ it[1] },
        ch_gtf.map{ it[1] },
        ch_fasta_fai.map{ it[1] },
        false
    )
    ch_seg_resolved_gtf = RESOLVE_UNANNOTATED.out.gtf
    ch_versions         = ch_versions.mix(RESOLVE_UNANNOTATED.out.versions)
    }

    /*
    * MODULE: Resolve the GTF regions that iCount did not annotate REGIONS FILE
    */
    if (params.regions_resolved_gtf){
    ch_regions_resolved_gtf = params.regions_resolved_gtf
    } else {
    RESOLVE_UNANNOTATED_REGIONS (
        ICOUNT_SEG_GTF.out.regions.map{ it[1] },
        ICOUNT_SEG_FILTGTF.out.regions.map{ it[1] },
        ch_gtf.map{ it[1] },
        ch_fasta_fai.map{ it[1] },
        false
    )
    ch_regions_resolved_gtf = RESOLVE_UNANNOTATED_REGIONS.out.gtf
    }

    /*
    * MODULE: Resolve the GTF regions that iCount did not annotate with genic_other flag
    */
    if (params.seg_resolved_gtf_genic){
    ch_seg_resolved_gtf_genic = params.seg_resolved_gtf_genic
    } else {
    RESOLVE_UNANNOTATED_GENIC_OTHER (
        ICOUNT_SEG_GTF.out.gtf.map{ it[1] },
        ICOUNT_SEG_FILTGTF.out.gtf.map{ it[1] },
        ch_gtf.map{ it[1] },
        ch_fasta_fai.map{ it[1] },
        true
    )
    ch_seg_resolved_gtf_genic = RESOLVE_UNANNOTATED_GENIC_OTHER.out.gtf
    }

    /*
    * MODULE: Resolve the GTF regions that iCount did not annotate with genic_other flag REGIONS FILE
    */
    if (params.regions_resolved_gtf_genic){
    ch_regions_resolved_gtf_genic = params.regions_resolved_gtf_genic
    } else {
    RESOLVE_UNANNOTATED_GENIC_OTHER_REGIONS (
        ICOUNT_SEG_GTF.out.regions.map{ it[1] },
        ICOUNT_SEG_FILTGTF.out.regions.map{ it[1] },
        ch_gtf.map{ it[1] },
        ch_fasta_fai.map{ it[1] },
        true
    )
    ch_regions_resolved_gtf_genic = RESOLVE_UNANNOTATED_GENIC_OTHER_REGIONS.out.gtf
    }


    emit:
    fasta                      = ch_fasta                  // channel: [ val(meta), [ fasta ] ]
    fasta_fai                  = ch_fasta_fai              // channel: [ val(meta), [ fai ] ]
    gtf                        = ch_gtf                    // channel: [ val(meta), [ gtf ] ]
    filtered_gtf               = ch_filt_gtf               // channel: [ val(meta), [ gtf ] ]
    chrom_sizes                = ch_chrom_sizes            // channel: [ val(meta), [ txt ] ]
    smrna_fasta                = ch_smrna_fasta            // channel: [ val(meta), [ fasta ] ]
    smrna_fasta_fai            = ch_smrna_fasta_fai        // channel: [ val(meta), [ fai ] ]
    smrna_chrom_sizes          = ch_smrna_chrom_sizes      // channel: [ val(meta), [ txt ] ]
    longest_transcript         = ch_longest_transcript     // channel: [ val(meta), [ txt ] ]
    longest_transcript_fai     = ch_longest_transcript_fai  // channel: [ val(meta), [ fai ] ]
    seg_gtf                    = ch_seg_gtf                // channel: [ val(meta), [ gtf ] ]
    seg_filt_gtf               = ch_seg_filt_gtf           // channel: [ val(meta), [ gtf ] ]
    seg_resolved_gtf           = ch_seg_resolved_gtf       // channel: [ val(meta), [ gtf ] ]
    seg_resolved_gtf_genic     = ch_seg_resolved_gtf_genic // channel: [ val(meta), [ gtf ] ]
    regions_gtf                = ch_regions_gtf                // channel: [ val(meta), [ gtf ] ]
    regions_filt_gtf           = ch_regions_filt_gtf           // channel: [ val(meta), [ gtf ] ]
    regions_resolved_gtf       = ch_regions_resolved_gtf       // channel: [ val(meta), [ gtf ] ]
    regions_resolved_gtf_genic = ch_regions_resolved_gtf_genic // channel: [ val(meta), [ gtf ] ]
    genome_index               = ch_genome_index           // channel: [ val(meta), [ star_index ] ]
    smrna_index                = ch_smrna_index            // channel: [ val(meta), [ bt2_index ] ]
    versions                   = ch_versions               // channel: [ versions.yml ]
}
