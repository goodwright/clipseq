//
// Prepare all genome files for running the clipseq analysis pipeline
//

/*
* MODULES
*/
include { GUNZIP                                        } from '../../../../modules/nf-core/gunzip/main.nf'
include { LINUX as ONE_LINE_FASTA                       } from '../../../../modules/goodwright/linux/command/main.nf'

/*
* SUBWORKFLOWS
*/
include { PREPARE_REF     as PREPARE_RRNA_SEQUENCE                   } from '../prepare_ref/main.nf'
include { PREPARE_REF     as PREPARE_CANONICAL_SNRNA_SEQUENCE        } from '../prepare_ref/main.nf'
include { PREPARE_REF     as PREPARE_MATURE_TRNA_SEQUENCE            } from '../prepare_ref/main.nf'
include { PREPARE_REF     as PREPARE_MATURE_SNRNA_SEQUENCE           } from '../prepare_ref/main.nf'
include { PREPARE_ALINGER as PREPARE_RRNA                            } from '../prepare_aligner/main.nf'
include { PREPARE_ALINGER as PREPARE_MATURE_TRNA                     } from '../prepare_aligner/main.nf'
include { PREPARE_ALINGER as PREPARE_IMMATURE_TRNA                   } from '../prepare_aligner/main.nf'
include { PREPARE_ALINGER as PREPARE_MATURE_SNRNA                    } from '../prepare_aligner/main.nf'
include { PREPARE_ALINGER as PREPARE_IMMATURE_SNRNA                  } from '../prepare_aligner/main.nf'
include { PREPARE_ALINGER as PREPARE_CANONICAL_SNRNA                 } from '../prepare_aligner/main.nf'
include { PREPARE_ALINGER as PREPARE_SCA_SNO_Y                       } from '../prepare_aligner/main.nf'
include { PREPARE_ALINGER as PREPARE_MITO                            } from '../prepare_aligner/main.nf'
include { PREPARE_ALINGER as PREPARE_REPEATS                         } from '../prepare_aligner/main.nf'
include { PREPARE_ALINGER as PREPARE_GENOME_MINUS_MITO               } from '../prepare_aligner/main.nf'



workflow PREPARE_NCRNA {
    take:
    species                          // channel: [ "species 2 letter code" ]
    fasta                            // channel: [ fasta ]
    gtf                              // channel: [ gtf ]
    rDNA_sequence                    // channel: [ fasta ]
    canonical_snRNA_sequences        // channel: [ fasta ]
    mature_tRNA_sequence             // channel: [ fasta ]
    mature_snRNA_sequence            // channel: [ fasta ]
    immature_tRNA_bed                // channel: [ tar.gz ]
    repeats                          // channel: [ txt.gz ]
    mito_chromosome                  // channel: [ "name of mitochondrial chromosome" ]
    rRNA_index                       // channel: [ folder/tar.gz ]
    mature_tRNA_index                // channel: [ folder/tar.gz ]
    immature_tRNA_index              // channel: [ folder/tar.gz ]
    mature_snRNA_index               // channel: [ folder/tar.gz ]
    immature_snRNA_index             // channel: [ folder/tar.gz ]
    canonical_snRNA_index            // channel: [ folder/tar.gz ]
    mito_index                       // channel: [ folder/tar.gz ]
    repeats_index                    // channel: [ folder/tar.gz ]
    genome_minus_mito_index          // channel: [ folder/tar.gz ]
    
    main:
    ch_versions = Channel.empty()



    /*
    * SUBWORKFLOW: Uncompress and prepare main sequence files
    */
    if (!rDNA_sequence.toString().endsWith(".gz")){
        ch_rDNA_sequence       = rDNA_sequence
    } else {
        PREPARE_RRNA_SEQUENCE (
            rDNA_sequence,
            [],
            [],
            []
        )
        ch_rDNA_sequence       = PREPARE_RRNA_SEQUENCE.out.fasta
        ch_versions    = ch_versions.mix(PREPARE_RRNA_SEQUENCE.out.versions)
    }

    if (!canonical_snRNA_sequences.toString().endsWith(".gz")){
        ch_canonical_snRNA_sequences      = canonical_snRNA_sequences
    } else {
        PREPARE_CANONICAL_SNRNA_SEQUENCE (
            canonical_snRNA_sequences,
            [],
            [],
            []
        )
        ch_canonical_snRNA_sequences       = PREPARE_CANONICAL_SNRNA_SEQUENCE.out.fasta
    }

    if (!mature_tRNA_sequence.toString().endsWith(".gz")){
        ch_mature_tRNA_sequence      = mature_tRNA_sequence
    } else {
        PREPARE_MATURE_TRNA_SEQUENCE (
            mature_tRNA_sequence,
            [],
            [],
            []
        )
        ch_mature_tRNA_sequence       = PREPARE_MATURE_TRNA_SEQUENCE.out.fasta
    }

    if (!mature_snRNA_sequence.toString().endsWith(".gz")){
        ch_mature_snRNA_sequence      = mature_snRNA_sequence
    } else {
        PREPARE_MATURE_SNRNA_SEQUENCE (
            mature_snRNA_sequence,
            [],
            [],
            []
        )
        ch_mature_snRNA_sequence       = PREPARE_MATURE_SNRNA_SEQUENCE.out.fasta
    }

    /*
    * SUBWORKFLOW: Prepare rRNA BT1 index
    */
    PREPARE_RRNA (
        ["bowtie"],
        ch_rDNA_sequence,
        [],
        rRNA_index,
        []
    )
    ch_rRNA_index  = PREPARE_RRNA.out.bt2_index
    ch_versions    = ch_versions.mix(PREPARE_RRNA.out.versions)

    /*
    * SUBWORKFLOW: Prepare canonical snRNA BT1 index
    */
    PREPARE_CANONICAL_SNRNA (
        ["bowtie"],
        ch_canonical_snRNA_sequences,
        [],
        canonical_snRNA_index,
        []
    )
    ch_canonical_snRNA_index = PREPARE_RRNA.out.bt2_index

    /*
    * Prepare all mature snRNA - we have to grep it out of the transcripts fasta file
    */
    ch_mature_snRNA_sequence = mature_snRNA_sequence
    if ( mature_snRNA_sequence.toString().endsWith("gz")){
        GUNZIP (
            ch_mature_snRNA_sequence
        )
        ch_mature_snRNA_sequence = GUNZIP.out.gunzip
        ch_versions              = ch_versions.mix(GUNZIP.out.versions)
    }

    ONE_LINE_FASTA (
        ch_mature_snRNA_sequence
    )
    ch_mature_snRNA_sequence = ONE_LINE_FASTA.out.file
    ch_versions              = ch_versions.mix(ONE_LINE_FASTA.out.versions)

    RETRIEVE_MATURE_SNRNA (
        ONE_LINE_FASTA.out.file
    )
    ch_mature_snRNA_sequence = RETRIEVE_MATURE_SNRNA.out.file
    ch_versions              = ch_versions.mix(RETRIEVE_MATURE_SNRNA.out.versions)

    PREPARE_MATURE_SNRNA (
        ["bowtie"],
        ch_mature_snRNA_sequence,
        [],
        mature_snRNA_index,
        []
    )
    ch_mature_snRNA_index = PREPARE_MATURE_SNRNA.out.bt2_index

    /*
    * Prepare immature snRNA - we have to get gene intervals, extend with flanks and get sequence from main genome
    */





    /*
    * SUBWORKFLOW: Prepare STAR index for primary genome
    */
    PREPARE_PRIMARY_INDEX (
        ["star"],
        ch_fasta,
        ch_gtf_with_meta,
        [],
        genome_index_path
    )
    ch_genome_index = PREPARE_PRIMARY_INDEX.out.star_index
    ch_versions     = ch_versions.mix(PREPARE_PRIMARY_INDEX.out.versions)








    emit:
    rRNA_index                = ch_rRNA_index 
    mature_tRNA_index         = ch_mature_tRNA_index 
    immature_tRNA_index       = ch_immature_tRNA_index 
    mature_snRNA_index        = ch_mature_snRNA_index 
    immature_snRNA_index      = ch_immature_snRNA_index 
    canonical_snRNA_index     = ch_canonical_snRNA_index 
    mito_index                = ch_mito_index 
    repeats_index             = ch_repeats_index 
    genome_minus_mito_index   = ch_genome_minus_mito_index 
}
