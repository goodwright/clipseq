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

include { PREPARE_CLIPSEQ } from './subworkflows/goodwright/prepare_genome/prepare_clipseq'

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
    ch_target_genome_index = Channel.empty()
    ch_smrna_genome_index  = Channel.empty()

    // Prepare manditory params into channels 
    ch_samplesheet = file(params.samplesheet, checkIfExists: true)
    ch_fasta       = file(params.fasta, checkIfExists: true)
    ch_smrna_fasta = file(params.smrna_fasta, checkIfExists: true)
    ch_gtf         = file(params.gtf, checkIfExists: true)

    // Pepare non-manditory params into channels
    if(params.target_genome_index) { ch_target_genome_index = file(params.target_genome_index, checkIfExists: true) }
    if(params.smrna_genome_index) { ch_smrna_genome_index = file(params.smrna_genome_index, checkIfExists: true) }

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

    // fasta                  = ch_fasta                  // channel: [ val(meta), [ fasta ] ]
    // fasta_fai              = ch_fasta_fai              // channel: [ val(meta), [ fai ] ]
    // gtf                    = ch_gtf                    // channel: [ val(meta), [ gtf ] ]
    // filtered_gtf           = ch_filt_gtf               // channel: [ val(meta), [ gtf ] ]
    // chrom_sizes            = ch_chrom_sizes            // channel: [ val(meta), [ txt ] ]
    // smrna_fasta            = ch_smrna_fasta            // channel: [ val(meta), [ fasta ] ]
    // smrna_fasta_fai        = ch_smrna_fasta_fai        // channel: [ val(meta), [ fai ] ]
    // smrna_chrom_sizes      = ch_smrna_chrom_sizes      // channel: [ val(meta), [ txt ] ]
    // longest_transcript     = ch_longest_transcript     // channel: [ val(meta), [ txt ] ]
    // seg_gtf                = ch_seg_gtf                // channel: [ val(meta), [ gtf ] ]
    // seg_filt_gtf           = ch_seg_filt_gtf           // channel: [ val(meta), [ gtf ] ]
    // seg_resolved_gtf       = ch_seg_resolved_gtf       // channel: [ val(meta), [ gtf ] ]
    // seg_resolved_gtf_genic = ch_seg_resolved_gtf_genic // channel: [ val(meta), [ gtf ] ]
    // versions               = ch_versions               // channel: [ versions.yml ]
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