#!/usr/bin/env nextflow

include { PREPARE_CLIPSEQ } from '../../subworkflows/goodwright/prepare_genome/prepare_clipseq/main.nf'

workflow  {
    fasta         = file(params.fasta, checkIfExists: true)
    smrna_fasta   = file(params.smrna_fasta, checkIfExists: true)
    gtf           = file(params.gtf, checkIfExists: true)

    star_index                      = params.star_index ? file(params.star_index, checkIfExists: true) : null
    bowtie2_index                   = params.bowtie2_index ? file(params.bowtie2_index, checkIfExists: true) : null

    fasta_fai                  = []
    filtered_gtf               = []
    chrom_sizes                = []
    smrna_fasta_fai            = []
    smrna_chrom_sizes          = []
    longest_transcript         = []
    seg_gtf                    = []
    seg_filt_gtf               = []
    seg_resolved_gtf           = []
    seg_resolved_gtf_genic     = []
    regions_gtf                = []
    regions_filt_gtf           = []
    regions_resolved_gtf       = []
    regions_resolved_gtf_genic = []
    longest_transcript_fai     = []
    longest_transcript_gtf     = []

    if(params.fasta_fai) { fasta_fai = Channel.of([[:],file(params.fasta_fai, checkIfExists: true)]) } 
    if(params.filtered_gtf) { filtered_gtf = Channel.of([[:],file(params.filtered_gtf, checkIfExists: true)]) }
    if(params.chrom_sizes) { chrom_sizes = Channel.of([[:],file(params.chrom_sizes, checkIfExists: true)]) }
    if(params.smrna_fasta_fai) { smrna_fasta_fai = Channel.of([[:],file(params.smrna_fasta_fai, checkIfExists: true)]) }
    if(params.smrna_chrom_sizes) { smrna_chrom_sizes = Channel.of([[:],file(params.smrna_chrom_sizes, checkIfExists: true)]) }
    if(params.longest_transcript) { longest_transcript = Channel.of([[:],file(params.longest_transcript, checkIfExists: true)]) }
    if(params.longest_transcript_fai) { longest_transcript_fai = Channel.of([[:],file(params.longest_transcript_fai, checkIfExists: true)]) }
    if(params.longest_transcript_gtf) { longest_transcript_gtf = Channel.of([[:],file(params.longest_transcript_gtf, checkIfExists: true)]) }
    if(params.seg_gtf) { seg_gtf = Channel.of([[:],file(params.seg_gtf, checkIfExists: true)]) }
    if(params.seg_filt_gtf) { seg_filt_gtf = Channel.of([[:],file(params.seg_filt_gtf, checkIfExists: true)]) }
    if(params.seg_resolved_gtf) { seg_resolved_gtf = file(params.seg_resolved_gtf, checkIfExists: true) }
    if(params.seg_resolved_gtf_genic) { seg_resolved_gtf_genic= Channel.of([[:],file(params.seg_resolved_gtf_genic, checkIfExists: true)]) }
    if(params.regions_gtf) { regions_gtf = Channel.of([[:],file(params.regions_gtf, checkIfExists: true)]) }
    if(params.regions_filt_gtf) { regions_filt_gtf = Channel.of([[:],file(params.regions_filt_gtf, checkIfExists: true)]) }
    if(params.regions_resolved_gtf) { regions_resolved_gtf = file(params.regions_resolved_gtf, checkIfExists: true) }
    if(params.regions_resolved_gtf_genic) { regions_resolved_gtf_genic = Channel.of([[:],file(params.regions_resolved_gtf_genic, checkIfExists: true)]) }
   
    PREPARE_CLIPSEQ (
        fasta,
        smrna_fasta,
        gtf,
        star_index,
        bowtie2_index,
        fasta_fai,
        filtered_gtf,
        chrom_sizes,
        smrna_fasta_fai,
        smrna_chrom_sizes,
        longest_transcript,
        seg_gtf,
        seg_filt_gtf,
        seg_resolved_gtf,
        seg_resolved_gtf_genic,
        regions_gtf,
        regions_filt_gtf,
        regions_resolved_gtf,
        regions_resolved_gtf_genic,
        longest_transcript_fai,
        longest_transcript_gtf
    )
}
