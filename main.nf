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

include { summary_log     } from './modules/goodwright/util/logging/main'
include { multiqc_summary } from './modules/goodwright/util/logging/main'

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
    params.smrna_genome_index,
    params.fasta_fai,
    params.filtered_gtf,
    params.chrom_sizes,
    params.smrna_fasta_fai,
    params.smrna_chrom_sizes,
    params.longest_transcript,
    params.longest_transcript_fai,
    params.longest_transcript_gtf,
    params.seg_gtf,
    params.seg_filt_gtf,
    params.seg_resolved_gtf,
    params.seg_resolved_gtf_genic,
    params.regions_gtf,
    params.regions_filt_gtf,
    params.regions_resolved_gtf,
    params.regions_resolved_gtf_genic
]
for (param in check_param_list) { if (param) { file(param, checkIfExists: true) } }


// Stage dummy file to be used as an optional input where required
ch_dummy_file = file("$projectDir/assets/dummy_file.txt", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULEs
//

include { MULTIQC } from './modules/local/multiqc'
include { GET_CROSSLINKS as CALC_SMRNA_K1_CROSSLINKS      } from './modules/local/get_crosslinks'
include { GET_CROSSLINKS as CALC_GENOME_CROSSLINKS        } from './modules/local/get_crosslinks'
include { GET_CROSSLINKS as CALC_TRANSCRIPT_CROSSLINKS    } from './modules/local/get_crosslinks'


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

include { SAMTOOLS_SIMPLE_VIEW as FILTER_TRANSCRIPTS } from './modules/goodwright/samtools/simple_view/main'
include { CLIPPY as CLIPPY_GENOME                    } from './modules/goodwright/clippy/main'
include { CLIPPY as CLIPPY_TRANSCRIPT                } from './modules/goodwright/clippy/main'
include { PEKA                                       } from './modules/goodwright/peka/main'
include { DUMP_SOFTWARE_VERSIONS                     } from './modules/goodwright/dump_software_versions/main'
include { CLIPSEQ_CLIPQC                             } from './modules/goodwright/clipseq/clipqc/main'
include { ENCODE_MOVEUMI                             } from './modules/goodwright/clipseq/encode_moveumi/main'

//
// SUBWORKFLOWS
//

include { PREPARE_CLIPSEQ                                       } from './subworkflows/goodwright/prepare_genome/prepare_clipseq/main'
include { PARSE_FASTQ_INPUT                                     } from './subworkflows/goodwright/parse_fastq_input/main'
include { FASTQC_TRIMGALORE                                     } from './subworkflows/goodwright/fastqc_trimgalore/main'
include { RNA_ALIGN                                             } from './subworkflows/goodwright/rna_align/main'
include { BAM_DEDUP_SAMTOOLS_UMITOOLS as GENOME_UNIQUE_DEDUP    } from './subworkflows/goodwright/bam_dedup_samtools_umitools/main'
include { BAM_DEDUP_SAMTOOLS_UMITOOLS as GENOME_MULTI_DEDUP     } from './subworkflows/goodwright/bam_dedup_samtools_umitools/main'
include { BAM_DEDUP_SAMTOOLS_UMITOOLS as SMRNA_DEDUP            } from './subworkflows/goodwright/bam_dedup_samtools_umitools/main'
include { BAM_DEDUP_SAMTOOLS_UMITOOLS as SMRNA_K1_DEDUP         } from './subworkflows/goodwright/bam_dedup_samtools_umitools/main'
include { BAM_DEDUP_SAMTOOLS_UMITOOLS as TRANSCRIPT_DEDUP       } from './subworkflows/goodwright/bam_dedup_samtools_umitools/main'
include { PARACLU_ANALYSE as PARACLU_ANALYSE_GENOME             } from './subworkflows/goodwright/paraclu_analyse/main'
include { PARACLU_ANALYSE as PARACLU_ANALYSE_TRANSCRIPT         } from './subworkflows/goodwright/paraclu_analyse/main'
include { ICOUNT_ANALYSE                                        } from './subworkflows/goodwright/icount_analyse/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULEs
//

include { UMITOOLS_EXTRACT                                 } from './modules/nf-core/umitools/extract/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_FILT_TRANSCRIPT   } from './modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_FILT_TRANSCRIPT } from './modules/nf-core/samtools/index/main'

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

    // Prepare manditory params into channels
    ch_samplesheet = file(params.samplesheet, checkIfExists: true)
    ch_fasta       = file(params.fasta, checkIfExists: true)
    ch_smrna_fasta = file(params.smrna_fasta, checkIfExists: true)
    ch_gtf         = file(params.gtf, checkIfExists: true)


    if (params.run_genome_prep) {
        // Prepare non-manditory params into channels
        ch_target_genome_index        = []
        ch_smrna_genome_index         = []
        ch_fasta_fai                  = []
        ch_filtered_gtf               = []
        ch_chrom_sizes                = []
        ch_smrna_fasta_fai            = []
        ch_smrna_chrom_sizes          = []
        ch_longest_transcript         = []
        ch_seg_gtf                    = []
        ch_seg_filt_gtf               = []
        ch_seg_resolved_gtf           = []
        ch_seg_resolved_gtf_genic     = []
        ch_regions_gtf                = []
        ch_regions_filt_gtf           = []
        ch_regions_resolved_gtf       = []
        ch_regions_resolved_gtf_genic = []
        ch_longest_transcript_fai     = []
        ch_longest_transcript_gtf     = []

        if(params.target_genome_index) { ch_target_genome_index = file(params.target_genome_index, checkIfExists: true) }
        if(params.smrna_genome_index)  { ch_smrna_genome_index = file(params.smrna_genome_index, checkIfExists: true) }
        if(params.fasta_fai) { ch_fasta_fai = Channel.of([[:],file(params.fasta_fai, checkIfExists: true)]) }
        if(params.filtered_gtf) { ch_filtered_gtf = Channel.of([[:],file(params.filtered_gtf, checkIfExists: true)]) }
        if(params.chrom_sizes) { ch_chrom_sizes = Channel.of([[:],file(params.chrom_sizes, checkIfExists: true)]) }
        if(params.smrna_fasta_fai) { ch_smrna_fasta_fai = Channel.of([[:],file(params.smrna_fasta_fai, checkIfExists: true)]) }
        if(params.smrna_chrom_sizes) { ch_smrna_chrom_sizes = Channel.of([[:],file(params.smrna_chrom_sizes, checkIfExists: true)]) }
        if(params.longest_transcript) { ch_longest_transcript = Channel.of([[:],file(params.longest_transcript, checkIfExists: true)]) }
        if(params.longest_transcript_fai) { ch_longest_transcript_fai = Channel.of([[:],file(params.longest_transcript_fai, checkIfExists: true)]) }
        if(params.longest_transcript_gtf) { ch_longest_transcript_gtf = Channel.of([[:],file(params.longest_transcript_gtf, checkIfExists: true)]) }
        if(params.seg_gtf) { ch_seg_gtf = Channel.of([[:],file(params.seg_gtf, checkIfExists: true)]) }
        if(params.seg_filt_gtf) { ch_seg_filt_gtf = Channel.of([[:],file(params.seg_filt_gtf, checkIfExists: true)]) }
        if(params.seg_resolved_gtf) { ch_seg_resolved_gtf = Channel.of([[:],file(params.seg_resolved_gtf, checkIfExists: true)]) }
        if(params.seg_resolved_gtf_genic) { ch_seg_resolved_gtf_genic= Channel.of([[:],file(params.seg_resolved_gtf_genic, checkIfExists: true)]) }
        if(params.regions_gtf) { ch_regions_gtf = Channel.of([[:],file(params.regions_gtf, checkIfExists: true)]) }
        if(params.regions_filt_gtf) { ch_regions_filt_gtf = Channel.of([[:],file(params.regions_filt_gtf, checkIfExists: true)]) }
        if(params.regions_resolved_gtf) { ch_regions_resolved_gtf = Channel.of([[:],file(params.regions_resolved_gtf, checkIfExists: true)]) }
        if(params.regions_resolved_gtf_genic) { ch_regions_resolved_gtf_genic = Channel.of([[:],file(params.regions_resolved_gtf_genic, checkIfExists: true)]) }

        /*
        * SUBWORKFLOW: Prepare clipseq genome files
        */
        PREPARE_CLIPSEQ (
            ch_fasta,
            ch_smrna_fasta,
            ch_gtf,
            ch_target_genome_index,
            ch_smrna_genome_index,
            ch_fasta_fai,
            ch_filtered_gtf,
            ch_chrom_sizes,
            ch_smrna_fasta_fai,
            ch_smrna_chrom_sizes,
            ch_longest_transcript,
            ch_seg_gtf,
            ch_seg_filt_gtf,
            ch_seg_resolved_gtf,
            ch_seg_resolved_gtf_genic,
            ch_regions_gtf,
            ch_regions_filt_gtf,
            ch_regions_resolved_gtf,
            ch_regions_resolved_gtf_genic,
            ch_longest_transcript_fai,
            ch_longest_transcript_gtf
        )
        ch_versions                   = ch_versions.mix(PREPARE_CLIPSEQ.out.versions)
        ch_fasta                      = PREPARE_CLIPSEQ.out.fasta
        ch_fasta_fai                  = PREPARE_CLIPSEQ.out.fasta_fai
        ch_gtf                        = PREPARE_CLIPSEQ.out.gtf
        ch_filtered_gtf               = PREPARE_CLIPSEQ.out.filtered_gtf
        ch_chrom_sizes                = PREPARE_CLIPSEQ.out.chrom_sizes
        ch_smrna_fasta                = PREPARE_CLIPSEQ.out.smrna_fasta
        ch_smrna_fasta_fai            = PREPARE_CLIPSEQ.out.smrna_fasta_fai
        ch_smrna_chrom_sizes          = PREPARE_CLIPSEQ.out.smrna_chrom_sizes
        ch_longest_transcript         = PREPARE_CLIPSEQ.out.longest_transcript
        ch_longest_transcript_fai     = PREPARE_CLIPSEQ.out.longest_transcript_fai
        ch_longest_transcript_gtf     = PREPARE_CLIPSEQ.out.longest_transcript_gtf
        ch_seg_gtf                    = PREPARE_CLIPSEQ.out.seg_gtf
        ch_seg_filt_gtf               = PREPARE_CLIPSEQ.out.seg_filt_gtf
        ch_seg_resolved_gtf           = PREPARE_CLIPSEQ.out.seg_resolved_gtf
        ch_seg_resolved_gtf_genic     = PREPARE_CLIPSEQ.out.seg_resolved_gtf_genic
        ch_regions_gtf                = PREPARE_CLIPSEQ.out.regions_gtf
        ch_regions_filt_gtf           = PREPARE_CLIPSEQ.out.regions_filt_gtf
        ch_regions_resolved_gtf       = PREPARE_CLIPSEQ.out.regions_resolved_gtf
        ch_regions_resolved_gtf_genic = PREPARE_CLIPSEQ.out.regions_resolved_gtf_genic
        ch_target_genome_index        = PREPARE_CLIPSEQ.out.genome_index
        ch_smrna_genome_index         = PREPARE_CLIPSEQ.out.smrna_index
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
    if(params.encode_eclip){
        ENCODE_MOVEUMI (
            ch_fastq
        )
        ch_versions = ch_versions.mix(ENCODE_MOVEUMI.out.versions)
        ch_fastq    = ENCODE_MOVEUMI.out.reads
    }
    if(params.run_move_umi_to_header){
        UMITOOLS_EXTRACT (
            ch_fastq
        )
        ch_versions = ch_versions.mix(UMITOOLS_EXTRACT.out.versions)
        ch_fastq    = UMITOOLS_EXTRACT.out.reads
    }

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

    ch_genome_bam     = Channel.empty()
    ch_genome_bai     = Channel.empty()
    ch_transcript_bam = Channel.empty()
    ch_transcript_bai = Channel.empty()
    ch_bt_log         = Channel.empty()
    ch_star_log       = Channel.empty()
    if(params.run_alignment) {
        /*
        * SUBWORKFLOW: Run alignment to target and smrna genome. sort/index the output
        */

        RNA_ALIGN (
            ch_fastq,
            ch_smrna_genome_index,
            ch_target_genome_index,
            ch_filtered_gtf,
            ch_fasta
        )
        ch_versions           = ch_versions.mix(RNA_ALIGN.out.versions)
        ch_genome_unique_bam  = RNA_ALIGN.out.genome_unique_bam
        ch_genome_unique_bai  = RNA_ALIGN.out.genome_unique_bai
        ch_genome_multi_bam   = RNA_ALIGN.out.genome_multi_bam
        ch_genome_multi_bai   = RNA_ALIGN.out.genome_multi_bai
        ch_smrna_bam          = RNA_ALIGN.out.smrna_bam
        ch_smrna_bai          = RNA_ALIGN.out.smrna_bai
        ch_smrna_k1_bam       = RNA_ALIGN.out.smrna_k1_bam
        ch_smrna_k1_bai       = RNA_ALIGN.out.smrna_k1_bai
        ch_transcript_bam     = RNA_ALIGN.out.transcript_bam
        ch_transcript_bai     = RNA_ALIGN.out.transcript_bai
        ch_bt_log             = RNA_ALIGN.out.bt_log
        ch_star_log           = RNA_ALIGN.out.star_log_final
    }

    if(params.run_read_filter) {
        /*
        * CHANNEL: Combine bam and bai files on id
        */
        ch_transcript_bam_bai = ch_transcript_bam
            .map { row -> [row[0].id, row ].flatten()}
            .join ( ch_transcript_bai.map { row -> [row[0].id, row ].flatten()} )
            .map { row -> [row[1], row[2], row[4]] }

        /*
        * MODULE: Filter transcriptome bam on longest transcripts
        */
        FILTER_TRANSCRIPTS (
            ch_transcript_bam_bai,
            [],
            ch_longest_transcript.collect{ it[1] }
        )
        ch_versions = ch_versions.mix(FILTER_TRANSCRIPTS.out.versions)

        /*
        * SUBWORKFLOW: sort, index filtered bam
        */
        SAMTOOLS_SORT_FILT_TRANSCRIPT ( FILTER_TRANSCRIPTS.out.bam )
        SAMTOOLS_INDEX_FILT_TRANSCRIPT ( SAMTOOLS_SORT_FILT_TRANSCRIPT.out.bam )

        ch_versions       = ch_versions.mix(SAMTOOLS_SORT_FILT_TRANSCRIPT.out.versions)
        ch_versions       = ch_versions.mix(SAMTOOLS_INDEX_FILT_TRANSCRIPT.out.versions)
        ch_transcript_bam = SAMTOOLS_SORT_FILT_TRANSCRIPT.out.bam
        ch_transcript_bai = SAMTOOLS_INDEX_FILT_TRANSCRIPT.out.bai
    }

    ch_umi_log = Channel.empty()
    if(params.run_umi_dedup) {
        /*
        * CHANNEL: Combine bam and bai files on id
        */
        ch_genome_unique_bam_bai = ch_genome_unique_bam
            .map { row -> [row[0].id, row ].flatten()}
            .join ( ch_genome_unique_bai.map { row -> [row[0].id, row ].flatten()} )
            .map { row -> [row[1], row[2], row[4]] }

        ch_genome_multi_bam_bai = ch_genome_multi_bam
            .map { row -> [row[0].id, row ].flatten()}
            .join ( ch_genome_multi_bai.map { row -> [row[0].id, row ].flatten()} )
            .map { row -> [row[1], row[2], row[4]] }

        ch_smrna_bam_bai = ch_smrna_bam
            .map { row -> [row[0].id, row ].flatten()}
            .join ( ch_smrna_bai.map { row -> [row[0].id, row ].flatten()} )
            .map { row -> [row[1], row[2], row[4]] }

        ch_smrna_k1_bam_bai = ch_smrna_k1_bam
            .map { row -> [row[0].id, row ].flatten()}
            .join ( ch_smrna_k1_bai.map { row -> [row[0].id, row ].flatten()} )
            .map { row -> [row[1], row[2], row[4]] }

        ch_transcript_bam_bai = ch_transcript_bam
            .map { row -> [row[0].id, row ].flatten()}
            .join ( ch_transcript_bai.map { row -> [row[0].id, row ].flatten()} )
            .map { row -> [row[1], row[2], row[4]] }

        /*
        * SUBWORKFLOW: Run umi deduplication on genome-level alignments
        */
        GENOME_UNIQUE_DEDUP (
            ch_genome_unique_bam_bai
        )
        ch_versions   = ch_versions.mix(GENOME_UNIQUE_DEDUP.out.versions)
        ch_genome_bam = GENOME_UNIQUE_DEDUP.out.bam
        ch_genome_bai = GENOME_UNIQUE_DEDUP.out.bai
        ch_umi_log    = GENOME_UNIQUE_DEDUP.out.umi_log

        GENOME_MULTI_DEDUP (
            ch_genome_multi_bam_bai
        )
        ch_versions   = ch_versions.mix(GENOME_MULTI_DEDUP.out.versions)

        SMRNA_DEDUP (
            ch_smrna_bam_bai
        )
        ch_versions   = ch_versions.mix(SMRNA_DEDUP.out.versions)

        SMRNA_K1_DEDUP (
            ch_smrna_k1_bam_bai
        )
        ch_versions     = ch_versions.mix(SMRNA_K1_DEDUP.out.versions)
        ch_smrna_k1_bam = SMRNA_K1_DEDUP.out.bam
        ch_smrna_k1_bai = SMRNA_K1_DEDUP.out.bai
        ch_umi_log      = SMRNA_K1_DEDUP.out.umi_log

        /*
        * SUBWORKFLOW: Run umi deduplication on transcript-level alignments
        */
        TRANSCRIPT_DEDUP (
            ch_transcript_bam_bai
        )
        ch_versions       = ch_versions.mix(TRANSCRIPT_DEDUP.out.versions)
        ch_transcript_bam = TRANSCRIPT_DEDUP.out.bam
        ch_transcript_bai = TRANSCRIPT_DEDUP.out.bai
    } else {
        ch_genome_bam = ch_genome_unique_bam
        ch_genome_bai = ch_genome_unique_bai
    }

    ch_genome_crosslink_bed           = Channel.empty()
    ch_genome_crosslink_coverage      = Channel.empty()
    ch_genome_crosslink_coverage_norm = Channel.empty()
    ch_trans_crosslink_bed            = Channel.empty()
    ch_trans_crosslink_coverage       = Channel.empty()
    ch_trans_crosslink_coverage_norm  = Channel.empty()
    if(params.run_calc_crosslinks) {
        /*
        * SUBWORKFLOW: Run crosslink calculation for smRNA with -k 1
        */
        CALC_SMRNA_K1_CROSSLINKS (
            ch_smrna_k1_bam.join(ch_smrna_k1_bai),
            ch_smrna_fasta_fai.collect{ it[1] },
            params.crosslink_position
        )
        ch_versions                      = ch_versions.mix(CALC_SMRNA_K1_CROSSLINKS.out.versions)
        ch_smrna_crosslink_bed           = CALC_SMRNA_K1_CROSSLINKS.out.bed
        ch_smrna_crosslink_coverage      = CALC_SMRNA_K1_CROSSLINKS.out.coverage
        ch_smrna_crosslink_coverage_norm = CALC_SMRNA_K1_CROSSLINKS.out.coverage_norm

        /*
        * SUBWORKFLOW: Run crosslink calculation for genome
        */
        CALC_GENOME_CROSSLINKS (
            ch_genome_bam.join(ch_genome_bai),
            ch_fasta_fai.collect{ it[1] },
            params.crosslink_position
        )
        ch_versions                       = ch_versions.mix(CALC_GENOME_CROSSLINKS.out.versions)
        ch_genome_crosslink_bed           = CALC_GENOME_CROSSLINKS.out.bed
        ch_genome_crosslink_coverage      = CALC_GENOME_CROSSLINKS.out.coverage
        ch_genome_crosslink_coverage_norm = CALC_GENOME_CROSSLINKS.out.coverage_norm

        /*
        * SUBWORKFLOW: Run crosslink calculation for transcripts
        */
        CALC_TRANSCRIPT_CROSSLINKS (
            ch_transcript_bam.join(ch_transcript_bai),
            ch_longest_transcript_fai.collect{ it[1] },
            params.crosslink_position
        )
        ch_versions                      = ch_versions.mix(CALC_TRANSCRIPT_CROSSLINKS.out.versions)
        ch_trans_crosslink_bed           = CALC_TRANSCRIPT_CROSSLINKS.out.bed
        ch_trans_crosslink_coverage      = CALC_TRANSCRIPT_CROSSLINKS.out.coverage
        ch_trans_crosslink_coverage_norm = CALC_TRANSCRIPT_CROSSLINKS.out.coverage_norm
    }

    if(params.run_peak_calling) {
        /*
        * MODULE: Run clippy on genome-level crosslinks
        */
        CLIPPY_GENOME (
            ch_genome_crosslink_bed,
            ch_filtered_gtf.collect{ it[1] },
            ch_fasta_fai.collect{ it[1] }
        )
        ch_versions = ch_versions.mix(CLIPPY_GENOME.out.versions)

        /*
        * MODULE: Run clippy on transcript-level crosslinks
        */
        CLIPPY_TRANSCRIPT (
            ch_trans_crosslink_bed,
            ch_longest_transcript_gtf.collect{ it[1] },
            ch_longest_transcript_fai.collect{ it[1] }
        )
        ch_versions = ch_versions.mix(CLIPPY_TRANSCRIPT.out.versions)

        /*
        * SUBWORKFLOW: Run paraclu on genome-level and transcript-level crosslinks
        */
        PARACLU_ANALYSE_GENOME (
            ch_genome_crosslink_bed,
            params.paraclu_min_value
        )
        ch_versions = ch_versions.mix(PARACLU_ANALYSE_GENOME.out.versions)

        PARACLU_ANALYSE_TRANSCRIPT (
            ch_trans_crosslink_bed,
            params.paraclu_min_value
        )

        /*
        * SUBWORKFLOW: Run iCount on genome-level crosslinks
        */
        ICOUNT_ANALYSE (
            ch_smrna_crosslink_bed,
            ch_genome_crosslink_bed,
            ch_regions_resolved_gtf.collect{ it[1] },
            ch_seg_resolved_gtf.collect{ it[1] },
            true
        )
        ch_versions = ch_versions.mix(ICOUNT_ANALYSE.out.versions)

        /*
        * MODULE: Run peka on genome-level crosslinks
        */
        ch_peka_input = CLIPPY_GENOME.out.peaks
            .join(ch_genome_crosslink_bed, by: 0)
        PEKA (
            ch_peka_input.map { meta, peaks, crosslinks -> [meta, peaks] },
            ch_peka_input.map { meta, peaks, crosslinks -> [meta, crosslinks] },
            ch_fasta.collect{ it[1] },
            ch_fasta_fai.collect{ it[1] },
            ch_regions_resolved_gtf.collect{ it[1] }
        )
        ch_versions = ch_versions.mix(PEKA.out.versions)
    }

    if(params.run_reporting) {
        /*
        * MODULE: Collect software versions
        */
        DUMP_SOFTWARE_VERSIONS (
            ch_versions.unique().collectFile()
        )

        /*
        * MODULE: Run clipqc
        */
        CLIPSEQ_CLIPQC (
            ch_bt_log.collect{ it[1] },
            ch_star_log.collect{ it[1] },
            ch_umi_log.collect{ it[1] },
            ch_genome_crosslink_bed.collect{ it[1] },
            ICOUNT_ANALYSE.out.bed_peaks.collect{ it[1] },
            PARACLU_ANALYSE_GENOME.out.peaks.collect{ it[1] },
            CLIPPY_GENOME.out.peaks.collect{ it[1] }
        )

        /*
        * MODULE: Run multiqc
        */
        workflow_summary    = multiqc_summary(workflow, params)
        ch_workflow_summary = Channel.value(workflow_summary)

        MULTIQC (
            ch_multiqc_config,
            DUMP_SOFTWARE_VERSIONS.out.mqc_yml.collect(),
            DUMP_SOFTWARE_VERSIONS.out.mqc_unique_yml.collect(),
            ch_workflow_summary.collectFile(name: "workflow_summary_mqc.yml"),
            FASTQC_TRIMGALORE.out.fastqc_zip.collect{it[1]}.ifEmpty([]),
            FASTQC_TRIMGALORE.out.fastqc_trim_zip.collect{it[1]}.ifEmpty([]),
            FASTQC_TRIMGALORE.out.trim_log.collect{it[1]}.ifEmpty([]),
            ch_bt_log.collect{it[1]}.ifEmpty([]),
            ch_star_log.collect{it[1]}.ifEmpty([]),
            CLIPSEQ_CLIPQC.out.tsv.collect().ifEmpty([])
        )
    }
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
