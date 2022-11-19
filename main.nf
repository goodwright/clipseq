#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    goodwright/clipseq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/goodwright/clipseq
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { summary_log } from './modules/goodwright/util/logging/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INIT
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

log.info summary_log(workflow, params, params.debug, params.monochrome_logs)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Check manditory input parameters to see if the files exist if they have been specified
check_param_list = [
    samplesheet: params.samplesheet,
    fasta: params.fasta,
    smrna_fasta: params.smrna_fasta,
    gtf: params.gtf
]
for (param in check_param_list) { 
    if (!param.value) { 
        exit 1, "Required parameter not specified: ${param.key}"
    } 
    else {
        file(param.value, checkIfExists: true)
    }
}

// Check non-manditory input parameters to see if the files exist if they have been specified
check_param_list = [
    params.target_genome_index,
    params.smrna_genome_index
]
for (param in check_param_list) { if (param) { file(param, checkIfExists: true) } }

// Stage dummy file to be used as an optional input where required
ch_dummy_file = file("$projectDir/assets/dummy_file.txt", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)

// // Header files for MultiQC
// ch_pca_header_multiqc        = file("$projectDir/assets/multiqc/deseq2_pca_header.txt", checkIfExists: true)
// ch_clustering_header_multiqc = file("$projectDir/assets/multiqc/deseq2_clustering_header.txt", checkIfExists: true)
// ch_biotypes_header_multiqc   = file("$projectDir/assets/multiqc/biotypes_header.txt", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULEs
//

//
// SUBWORKFLOWS
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT GOODWRIGHT MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULEs
//

//
// SUBWORKFLOWS
//

include { PREPARE_CLIPSEQ   } from './subworkflows/goodwright/prepare_genome/prepare_clipseq'
include { PARSE_FASTQ_INPUT } from './subworkflows/goodwright/parse_fastq_input'
include { FASTQC_TRIMGALORE } from './subworkflows/goodwright/fastqc_trimgalore/main'
include { RNA_ALIGN         } from './subworkflows/goodwright/rna_align/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULEs
//

//
// SUBWORKFLOWS
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CLIPSEQ {
    // Init
    ch_versions            = Channel.empty()
    ch_target_genome_index = []
    ch_smrna_genome_index  = []

    // Prepare manditory params into channels 
    ch_samplesheet = file(params.samplesheet, checkIfExists: true)
    ch_fasta       = file(params.fasta, checkIfExists: true)
    ch_smrna_fasta = file(params.smrna_fasta, checkIfExists: true)
    ch_gtf         = file(params.gtf, checkIfExists: true)

    // Pepare non-manditory params into channels
    if(params.target_genome_index) { ch_target_genome_index = file(params.target_genome_index, checkIfExists: true) }
    if(params.smrna_genome_index)  { ch_smrna_genome_index = file(params.smrna_genome_index, checkIfExists: true) }

    // Prepare genome and build indexes if required
    ch_fasta_fai              = Channel.empty()
    ch_filtered_gtf           = Channel.empty()
    ch_chrom_sizes            = Channel.empty()
    ch_smrna_fasta_fai        = Channel.empty()
    ch_smrna_chrom_sizes      = Channel.empty()
    ch_longest_transcript     = Channel.empty()
    ch_seg_gtf                = Channel.empty()
    ch_seg_filt_gtf           = Channel.empty()
    ch_seg_resolved_gtf       = Channel.empty()
    ch_seg_resolved_gtf_genic = Channel.empty()
    if (params.run_genome_prep) {
        /*
        * SUBWORKFLOW: Prepare clipseq genome files
        */
        PREPARE_CLIPSEQ (
            ch_fasta,
            ch_smrna_fasta,
            ch_gtf,
            ch_target_genome_index,
            ch_smrna_genome_index
        )
        ch_versions               = ch_versions.mix(PREPARE_CLIPSEQ.out.versions)
        ch_fasta                  = PREPARE_CLIPSEQ.out.fasta
        ch_fasta_fai              = PREPARE_CLIPSEQ.out.fasta_fai
        ch_gtf                    = PREPARE_CLIPSEQ.out.gtf
        ch_filtered_gtf           = PREPARE_CLIPSEQ.out.filtered_gtf
        ch_chrom_sizes            = PREPARE_CLIPSEQ.out.chrom_sizes
        ch_smrna_fasta            = PREPARE_CLIPSEQ.out.smrna_fasta
        ch_smrna_fasta_fai        = PREPARE_CLIPSEQ.out.smrna_fasta_fai
        ch_smrna_chrom_sizes      = PREPARE_CLIPSEQ.out.smrna_chrom_sizes
        ch_longest_transcript     = PREPARE_CLIPSEQ.out.longest_transcript
        ch_seg_gtf                = PREPARE_CLIPSEQ.out.seg_gtf
        ch_seg_filt_gtf           = PREPARE_CLIPSEQ.out.seg_filt_gtf
        ch_seg_resolved_gtf       = PREPARE_CLIPSEQ.out.seg_resolved_gtf
        ch_seg_resolved_gtf_genic = PREPARE_CLIPSEQ.out.seg_resolved_gtf_genic
        ch_target_genome_index    = PREPARE_CLIPSEQ.out.genome_index
        ch_smrna_genome_index     = PREPARE_CLIPSEQ.out.smrna_index
    }

    ch_fastq = Channel.empty()
    if(params.run_input_check) {
        /*
        * SUBWORKFLOW: Read in samplesheet, validate, stage input files and merge replicates
        */
        PARSE_FASTQ_INPUT (
            ch_samplesheet
        )
        ch_versions = ch_versions.mix(PARSE_FASTQ_INPUT.out.versions)
        ch_fastq    = PARSE_FASTQ_INPUT.out.fastq
    }
    //EXAMPLE CHANNEL STRUCT: [[id:h3k27me3_R1, group:h3k27me3, replicate:1, single_end:false], [FASTQ]]
    //ch_fastq | view

    if(params.run_trim_galore_fastqc) {
        /*
        * SUBWORKFLOW: Run fastqc and trimming
        */
        FASTQC_TRIMGALORE (
            ch_fastq,
            params.skip_fastqc,
            params.skip_trimming
        )
        ch_versions = ch_versions.mix(FASTQC_TRIMGALORE.out.versions)
        ch_fastq    = FASTQC_TRIMGALORE.out.fastq
    }
    //EXAMPLE CHANNEL STRUCT: [[id:h3k27me3_R1, group:h3k27me3, replicate:1, single_end:false], [FASTQ]]
    //ch_fastq | view

    if(params.run_alignment) {
        /*
        * SUBWORKFLOW: Run alignment to target and smrna genome. sort/index/stats the output
        */
        RNA_ALIGN (
            ch_fastq,
            ch_smrna_genome_index.map{ [ [id:'smrna' ], it[1] ] },
            ch_target_genome_index,
            ch_filtered_gtf.map{ it[1] },
            ch_fasta.map{ it[1] }
        )
        ch_versions                     = ch_versions.mix(RNA_ALIGN.out.versions)
        ch_genome_bam                   = RNA_ALIGN.out.genome_bam
        ch_genome_bai                   = RNA_ALIGN.out.genome_bai
        ch_genome_samtools_stats        = RNA_ALIGN.out.genome_stats
        ch_genome_samtools_flagstat     = RNA_ALIGN.out.genome_flagstat
        ch_genome_samtools_idxstats     = RNA_ALIGN.out.genome_idxstats
        ch_transcript_bam               = RNA_ALIGN.out.transcript_bam
        ch_transcript_bai               = RNA_ALIGN.out.transcript_bai
        ch_transcript_samtools_stats    = RNA_ALIGN.out.transcript_stats
        ch_transcript_samtools_flagstat = RNA_ALIGN.out.transcript_flagstat
        ch_transcript_samtools_idxstats = RNA_ALIGN.out.transcript_idxstats
    }


    // TODO: clipseq qc
    // TODO: software versions
    // TODO: multiqc

    // TODO: test list
    //     - tar indexes
    //     - replciates to merge
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    CLIPSEQ ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EVENTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// workflow.onComplete {
// }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/