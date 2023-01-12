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
include { CLIPPY                                     } from './modules/goodwright/clippy/main'
include { PEKA                                       } from './modules/goodwright/peka/main'
include { DUMP_SOFTWARE_VERSIONS                     } from './modules/goodwright/dump_software_versions/main'
include { CLIPSEQ_CLIPQC                             } from './modules/goodwright/clipseq/clipqc/main'

//
// SUBWORKFLOWS
//

include { PREPARE_CLIPSEQ                                       } from './subworkflows/goodwright/prepare_genome/prepare_clipseq/main'
include { PARSE_FASTQ_INPUT                                     } from './subworkflows/goodwright/parse_fastq_input/main'
include { FASTQC_TRIMGALORE                                     } from './subworkflows/goodwright/fastqc_trimgalore/main'
include { RNA_ALIGN                                             } from './subworkflows/goodwright/rna_align/main'
include { BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS as GENOME_DEDUP     } from './subworkflows/goodwright/bam_dedup_stats_samtools_umitools/main'
include { BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS as TRANSCRIPT_DEDUP } from './subworkflows/goodwright/bam_dedup_stats_samtools_umitools/main'
include { CLIP_CALC_CROSSLINKS as CALC_GENOME_CROSSLINKS        } from './subworkflows/goodwright/clip_calc_crosslinks/main'
include { CLIP_CALC_CROSSLINKS as CALC_TRANSCRIPT_CROSSLINKS    } from './subworkflows/goodwright/clip_calc_crosslinks/main'
include { PARACLU_ANALYSE                                       } from './subworkflows/goodwright/paraclu_analyse/main'
include { ICOUNT_ANALYSE                                        } from './subworkflows/goodwright/icount_analyse/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULEs
//

include { UMITOOLS_EXTRACT                           } from './modules/nf-core/umitools/extract/main'

//
// SUBWORKFLOWS
//

include { BAM_SORT_STATS_SAMTOOLS as BAM_SORT_STATS_SAMTOOLS_TRANSCRIPT } from './subworkflows/nf-core/bam_sort_stats_samtools/main'

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

    // Prepare non-manditory params into channels
    if(params.target_genome_index) { ch_target_genome_index = file(params.target_genome_index, checkIfExists: true) }
    if(params.smrna_genome_index)  { ch_smrna_genome_index = file(params.smrna_genome_index, checkIfExists: true) }

    // Prepare genome and build indexes if required
    ch_fasta_fai                  = Channel.empty()
    ch_filtered_gtf               = Channel.empty()
    ch_chrom_sizes                = Channel.empty()
    ch_smrna_fasta_fai            = Channel.empty()
    ch_smrna_chrom_sizes          = Channel.empty()
    ch_longest_transcript         = Channel.empty()
    ch_seg_gtf                    = Channel.empty()
    ch_seg_filt_gtf               = Channel.empty()
    ch_seg_resolved_gtf           = Channel.empty()
    ch_seg_resolved_gtf_genic     = Channel.empty()
	ch_regions_gtf                = Channel.empty()
    ch_regions_filt_gtf           = Channel.empty()
    ch_regions_resolved_gtf       = Channel.empty()
    ch_regions_resolved_gtf_genic = Channel.empty()

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

    ch_genome_bam                   = Channel.empty()
    ch_genome_bai                   = Channel.empty()
    ch_genome_samtools_stats        = Channel.empty()
    ch_genome_samtools_flagstat     = Channel.empty()
    ch_genome_samtools_idxstats     = Channel.empty()
    ch_transcript_bam               = Channel.empty()
    ch_transcript_bai               = Channel.empty()
    ch_transcript_samtools_stats    = Channel.empty()
    ch_transcript_samtools_flagstat = Channel.empty()
    ch_transcript_samtools_idxstats = Channel.empty()
    ch_bt_log                       = Channel.empty()
    ch_star_log                     = Channel.empty()
    if(params.run_alignment) {
        /*
        * SUBWORKFLOW: Run alignment to target and smrna genome. sort/index/stats the output
        */
        RNA_ALIGN (
            ch_fastq,
            ch_smrna_genome_index,
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
        ch_bt_log                       = RNA_ALIGN.out.bt_log
        ch_star_log                     = RNA_ALIGN.out.star_log_final
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
        * SUBWORKFLOW: Sort, index, stats on filtered bam
        */
        BAM_SORT_STATS_SAMTOOLS_TRANSCRIPT (
            FILTER_TRANSCRIPTS.out.bam,
            ch_fasta.map{ it[1] }
        )
        ch_versions                     = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS_TRANSCRIPT.out.versions)
        ch_transcript_bam               = BAM_SORT_STATS_SAMTOOLS_TRANSCRIPT.out.bam
        ch_transcript_bai               = BAM_SORT_STATS_SAMTOOLS_TRANSCRIPT.out.bai
        ch_transcript_samtools_stats    = BAM_SORT_STATS_SAMTOOLS_TRANSCRIPT.out.stats
        ch_transcript_samtools_flagstat = BAM_SORT_STATS_SAMTOOLS_TRANSCRIPT.out.flagstat
        ch_transcript_samtools_idxstats = BAM_SORT_STATS_SAMTOOLS_TRANSCRIPT.out.idxstats
    }

    ch_umi_log = Channel.empty()
    if(params.run_umi_dedup) {
        /*
        * CHANNEL: Combine bam and bai files on id
        */
        ch_genome_bam_bai = ch_genome_bam
            .map { row -> [row[0].id, row ].flatten()}
            .join ( ch_genome_bai.map { row -> [row[0].id, row ].flatten()} )
            .map { row -> [row[1], row[2], row[4]] }

        ch_transcript_bam_bai = ch_transcript_bam
            .map { row -> [row[0].id, row ].flatten()}
            .join ( ch_transcript_bai.map { row -> [row[0].id, row ].flatten()} )
            .map { row -> [row[1], row[2], row[4]] }

        /*
        * SUBWORKFLOW: Run umi deduplication on genome-level alignments
        */
        GENOME_DEDUP (
            ch_genome_bam_bai
        )
        ch_versions                 = ch_versions.mix(GENOME_DEDUP.out.versions)
        ch_genome_bam               = GENOME_DEDUP.out.bam
        ch_genome_bai               = GENOME_DEDUP.out.bai
        ch_genome_samtools_stats    = GENOME_DEDUP.out.stats
        ch_genome_samtools_flagstat = GENOME_DEDUP.out.flagstat
        ch_genome_samtools_idxstats = GENOME_DEDUP.out.idxstats
        ch_umi_log                  = GENOME_DEDUP.out.umi_log

        /*
        * SUBWORKFLOW: Run umi deduplication on transcript-level alignments
        */
        TRANSCRIPT_DEDUP (
            ch_transcript_bam_bai
        )
        ch_versions                     = ch_versions.mix(TRANSCRIPT_DEDUP.out.versions)
        ch_transcript_bam               = TRANSCRIPT_DEDUP.out.bam
        ch_transcript_bai               = TRANSCRIPT_DEDUP.out.bai
        ch_transcript_samtools_stats    = TRANSCRIPT_DEDUP.out.stats
        ch_transcript_samtools_flagstat = TRANSCRIPT_DEDUP.out.flagstat
        ch_transcript_samtools_idxstats = TRANSCRIPT_DEDUP.out.idxstats
    }

    ch_genome_crosslink_bed           = Channel.empty()
    ch_genome_crosslink_coverage      = Channel.empty()
    ch_genome_crosslink_coverage_norm = Channel.empty()
    ch_trans_crosslink_bed            = Channel.empty()
    ch_trans_crosslink_coverage       = Channel.empty()
    ch_trans_crosslink_coverage_norm  = Channel.empty()
    if(params.run_calc_crosslinks) {
        /*
        * SUBWORKFLOW: Run crosslink calculation for genome
        */
        CALC_GENOME_CROSSLINKS (
            ch_genome_bam,
            ch_fasta_fai.collect{ it[1] }
        )
        ch_versions                       = ch_versions.mix(CALC_GENOME_CROSSLINKS.out.versions)
        ch_genome_crosslink_bed           = CALC_GENOME_CROSSLINKS.out.bed
        ch_genome_crosslink_coverage      = CALC_GENOME_CROSSLINKS.out.coverage
        ch_genome_crosslink_coverage_norm = CALC_GENOME_CROSSLINKS.out.coverage_norm

        /*
        * SUBWORKFLOW: Run crosslink calculation for transcripts
        */
        CALC_TRANSCRIPT_CROSSLINKS (
            ch_genome_bam,
            ch_fasta_fai.collect{ it[1] }
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
        CLIPPY (
            ch_genome_crosslink_bed,
            ch_filtered_gtf.collect{ it[1] },
            ch_fasta_fai.collect{ it[1] }
        )
        ch_versions = ch_versions.mix(CLIPPY.out.versions)

        /*
        * SUBWORKFLOW: Run paraclu on genome-level crosslinks
        */
        PARACLU_ANALYSE (
            ch_genome_crosslink_bed,
            params.paraclu_min_value
        )
        ch_versions = ch_versions.mix(PARACLU_ANALYSE.out.versions)

        /*
        * SUBWORKFLOW: Run iCount on genome-level crosslinks
        */
        ICOUNT_ANALYSE (
            ch_genome_crosslink_bed,
            ch_seg_gtf.collect{ it[1] },
            ch_seg_resolved_gtf,
            true
        )
        ch_versions = ch_versions.mix(ICOUNT_ANALYSE.out.versions)

        // /*
        // * MODULE: Run peka on genome-level crosslinks
        // */
        PEKA (
            CLIPPY.out.peaks,
            ch_genome_crosslink_bed,
            ch_fasta.collect{ it[1] },
            ch_fasta_fai.collect{ it[1] },
            ch_regions_resolved_gtf
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
            ch_bt_log.map{ it[1] },
            ch_star_log.map{ it[1] },
            ch_umi_log.map{ it[1] },
            ch_genome_crosslink_bed.map{ it[1] },
            ICOUNT_ANALYSE.out.bed_peaks.map{ it[1] },
            PARACLU_ANALYSE.out.peaks.map{ it[1] },
            CLIPPY.out.peaks.map{ it[1] }
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
