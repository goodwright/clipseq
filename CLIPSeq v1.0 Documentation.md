# CLIPSeq v1.0 📎 

Full source code available on [Github](https://github.com/goodwright/clipseq). Specific instructions for Flow usage can be found [**here**]. 

## Pipeline Summary

Cross-linking and immunoprecipitation followed by sequencing (CLIP) has allowed high resolution studies of RNA binding protein (RBP)-RNA interactions at transcriptomic scale.

This pipeline enables analysis of various forms of single-end CLIP data including variants of iCLIP (eg. irCLIP, iCLIP2, iiCLIP) and eCLIP (note that we have achieved comparable results using this pipeline to study reformatted paired-end eCLIP also). The pipeline currently doesn't support mutation calling and therefore might not be suitable for PAR-CLIP analysis, but we plan to include this at a future date.

Reads are first 3' end trimmed for adapters and sequence quality before being mapped using Bowtie to a "small RNA" index consisting of rRNA and tRNA sequences. Reads that remain unmapped are then mapped to the chosen species' genome. Aligned reads are then deduplicated according to unique molecular identifier (UMI) and start position. Start position minus one of deduplicated reads is taken as the crosslinking position and crosslink bed files are produced. Peaks are then called using several tools: Paraclu, Clippy and iCount. Summary files reporting transcriptomic region, RNA biotype and gene level crosslink counts are produced. RNA maps of crosslinks are drawn around important transcriptomic landmarks. Crosslinks and Clippy peaks are taken as input to PEKA, to determine postionally enriched kmer motifs. Finally a summary of quality control metrics is generated and presented as a MultiQC report.

For a detailed review of the considerations behind each analysis step you can consult "[Data Science Issues in Studying Protein–RNA Interactions with CLIP Technologies](https://www.annualreviews.org/doi/abs/10.1146/annurev-biodatasci-080917-013525)".

## Required Input

This pipeline requires demultiplexed fastq sample input and an associated metadata spreadsheet. We recommend using [our demultiplex pipeline](https://github.com/goodwright/flow-nf/tree/master/subworkflows/goodwright/clip_demultiplex) to produce these individual samples fastq.

## Commonly used parameters

**Moving UMI from fastq reads to read header**

If you are analysing demultiplexed sample files, then depending on where you have sourced your fastq from, the UMI and experimental barcode might still be present at the 5' end of reads, which will cause errors in downstream analysis. Note this is common for historical public data downloaded from ArrayExpress, for example. 

**Changing the UMI delimiter**

Depending on how you demultiplexed your reads, you may need to change the delimiter that UMICollapse uses to search for UMIs in your aligned reads. We default to using "rbc:" which is the output of Ultraplex, the demultiplexer we use in our demultiplexing pipeline.

**Turning off deduplication**

For older data with very short or no UMIs at all, you may want to skip the UMI deduplication step.

**Changing peak calling or PEKA parameters** 

When you are working with data you're already familiar with you might have specific parameters in mind for peak calling or PEKA.  

## Pipeline in Detail

This detailed description will present each <u>subworkflow</u>/*module* run in the pipeline and give detailed information of inputs and outputs with filenames, and default parameters.

1. <u>PARSE_FASTQ_INPUT</u>

   Checks samplesheet input, then formats it and fastq files into an appropriate channel format for Nextflow.

   i. *SAMPLE_BASE_SAMPLESHEET_CHECK*: This module checks that the samplesheet input meets requirements. It outputs a sanitised samplesheet in a ".csv" file.

   ii. *CAT_FASTQ*: This module concatenates fastq files labelled as being the same sample.

   The output of this subworkflow is a Nextflow channel containing a tuple of metadata and fastq files.

2. (Optional) <u>UMI_TO_HEADER</u>

3. <u>FASTQC_TRIMGALORE</u>

   i. *FASTQC*: Checks sequencing quality and adapter content of raw reads. Input - .Output - .

   ii. *TRIMGALORE*

4. <u>RNA_ALIGN</u>
5. <u>FILTER_TRANSCRIPTS</u>
6. <u>BAM_SORT_STATS_SAMTOOLS_TRANSCRIPT</u>
7. <u>GENOME_DEDUP</u>
8. <u>TRANSCRIPT_DEDUP</u>
9. <u>CALC_GENOME_CROSSLINKS</u>
10. <u>CALC_TRANSCRIPT_CROSSLINKS</u>
11. *CLIPPY*
12. <u>PARACLU_ANALYSE</u>
13. <u>ICOUNT_ANALYSE</u>
14. *PEKA*
15. <u>CLIPSEQ_CLIPQC</u>
16. *MULTIQC*



## FAQ

**What is the difference between this pipeline and nf-core/clipseq?**

This pipeline is based on nf-core/clipseq but is updated in a few ways. For one, the Nextflow code has been updated from DSL1 to the modular DSL2. This pipeline also includes additional tools, namely [PEKA](https://github.com/ulelab/peka) for positional motif analysis and [Clippy](https://github.com/ulelab/clippy) for peak calling.

**I would like to contribute to this pipeline, or make a request for a new feature. Where can I do this?**

## Authors

This DSL2 Nextflow pipeline is maintained by Goodwright. It was updated from the DSL1 [nf-core/clipseq](https://nf-co.re/clipseq) pipeline in collaboration with the original authors and Prof. Jernej Ule. 

## References