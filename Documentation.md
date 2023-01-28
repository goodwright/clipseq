# CLIPSeq v1.0 📎

Full source code available on [Github](https://github.com/goodwright/clipseq). Specific instructions for Flow usage can be found [**here**].

## Pipeline Summary

Cross-linking and immunoprecipitation followed by sequencing (CLIP) has allowed high resolution studies of RNA binding protein (RBP)-RNA interactions at transcriptomic scale.

This pipeline enables analysis of various forms of single-end CLIP data including variants of iCLIP (eg. irCLIP, iCLIP2, iiCLIP) and eCLIP (note that we have achieved comparable results using this pipeline to study reformatted paired-end eCLIP also). The pipeline currently doesn't support mutation calling and therefore might not be suitable for PAR-CLIP analysis, but we plan to include this at a future date.

Reads are first 3' end trimmed for adapters and sequence quality before being mapped using Bowtie to a "small RNA" index consisting of rRNA and tRNA sequences. Reads that remain unmapped are then mapped to the chosen species' genome. Aligned reads are then deduplicated according to unique molecular identifier (UMI) and start position. Start position minus one of deduplicated reads is taken as the crosslinking position and crosslink bed files are produced. Peaks are then called using several tools: Paraclu, Clippy and iCount. Summary files reporting transcriptomic region, RNA biotype and gene level crosslink counts are produced. RNA maps of crosslinks are drawn around important transcriptomic landmarks. Crosslinks and Clippy peaks are taken as input to PEKA, to determine postionally enriched kmer motifs. Finally a summary of quality control metrics is generated and presented as a MultiQC report.

For a detailed review of the considerations behind each analysis step you can consult "[Data Science Issues in Studying Protein–RNA Interactions with CLIP Technologies](https://www.annualreviews.org/doi/abs/10.1146/annurev-biodatasci-080917-013525)".

## Required Input

This pipeline requires demultiplexed fastq sample input and an associated metadata spreadsheet. We recommend using [our demultiplex pipeline](https://github.com/goodwright/flow-nf/tree/master/subworkflows/goodwright/clip_demultiplex) to produce these individual samples fastq.

## Curated Outputs

These are what we would consider to be the most commonly used outputs of the pipeline, so on the Flow platform we present these as "curated ouputs". All of the outputs are listed in the "Pipeline in Detail" section below.

## **Crosslink files**

## **Peak files**

## **Summary/analysis of crosslinks and peaks**

## Commonly used parameters

**Moving UMI from fastq reads to read header**

If you are analysing demultiplexed sample files, then depending on where you have sourced your fastq from, the UMI and experimental barcode might still be present at the 5' end of reads, which will cause errors in downstream analysis. Note this is common for historical public data downloaded from ArrayExpress, for example.

To enable this option set `run_move_umi_to_header = true` and ensure you provide the UMI format to `move_umi`, for example `move_umi='NNNNNN'`.

**Changing the UMI delimiter**

Depending on how you demultiplexed your reads, you may need to change the delimiter that UMICollapse uses to search for UMIs in your aligned reads. We default to using "rbc:" which is the output of Ultraplex, the demultiplexer we use in our demultiplexing pipeline.

To change, set `umi_separator`.

**Turning off deduplication**

For older data with very short or no UMIs at all, you may want to skip the UMI deduplication step.

To turn of deduplications, set `run_umi_dedup = false`

**Changing other individual tool parameters**

When you are working with data you're already familiar with you might have specific parameters in mind for certain tools, here are your options for changing:

- Minimum length of reads kept by Trim Galore! after trimming, `trim_length`, eg. `10`
- Bowtie parameters for pre-mapping, `bowtie_params`, eg. `"-v 2 -m 100 --norc --best --strata"`
- Paraclu minimum cut off value, `paraclu_min_value`, eg. `30`

## Pipeline in Detail

This detailed description will present each <u>subworkflow</u>/_module_ run in the pipeline and give detailed information of inputs and outputs with filenames, and default parameters.

0. <u>PREPARE_CLIPSEQ</u>

   Prepares all of the alignment indexes and annotation files required by the pipeline, this only needs to be run once per genome, after this the output files can be provided to the pipeline to avoid this subworkflow running.

1. <u>PARSE_FASTQ_INPUT</u>

   Checks samplesheet input, then formats it and fastq files into an appropriate channel format for Nextflow.

   i. _SAMPLE_BASE_SAMPLESHEET_CHECK_: This module checks that the samplesheet input meets requirements. It outputs a sanitised samplesheet in a ".csv" file.

   ii. _CAT_FASTQ_: This module concatenates fastq files labelled as being the same sample.

   The output of this subworkflow is a Nextflow channel containing a tuple of metadata and fastq files.

2. (Optional) <u>UMI_TO_HEADER</u>

3. <u>FASTQC_TRIMGALORE</u>

   i. _FASTQC_: Checks sequencing quality and adapter content of raw reads.

   **⚪Input - sample `fastq`. 🟢Output - sample_fastqc.html `fastqc_html`, sample_fastqc.zip `fastqc_zip`**.

   ii. _TRIMGALORE_: Trims for 3'adapters, read quality at 3' end and removes reads that become too short after this. [Manual.](https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md)

   **🟣Parameters - `--fastqc --length ${params.trim_length} -q20`** Default params.trim_length = 10.

   **⚪Input - sample `fastq`. 🟢Output - sample_trimmed.fq.gz `fastq`, sample_trimming_report.txt `trim_log`, sample_trimmed_fastqc.html `fastqc_trim_html`, sample_trimmed_fastqc.zip `fastqc_trim_zip`**.

4. <u>RNA_ALIGN</u>

   i. _BOWTIE_ALIGN_ : Aligns reads to rRNA, tRNA sequences referred to as small RNA, smRNA. [Manual.](https://bowtie-bio.sourceforge.net/manual.shtml)

   **🟣Parameters - `${params.bowtie_params}`** Default params.bowtie_params = "-v 2 -m 100 --norc --best --strata"

   **⚪Input - sample_trimmed.fq.gz `reads`, bowtie (folder with bowtie index) `index`. 🟢Output - sample.bam `bam`, sample.out `log`, sample.unmapped.fastq.gz `fastq`.**

   ii. _STAR_ALIGN_ : Aligns unmapped reads from smRNA Bowtie mapping to the genome, also referred to as target. [Manual.](<[https://bowtie-bio.sourceforge.net/manual.shtml](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)>)

   **🟣Parameters - `--readFilesCommand zcat 
         --outSAMtype BAM SortedByCoordinate
         --quantMode TranscriptomeSAM ${params.star_params}`** Default params.star_params = "
   --outFilterMultimapNmax 1
   --outFilterMultimapScoreRange 1
   --outSAMattributes All
   --alignSJoverhangMin 8
   --alignSJDBoverhangMin 1
   --outFilterType BySJout
   --alignIntronMin 20
   --alignIntronMax 1000000
   --outFilterScoreMin 10
   --alignEndsType Extend5pOfRead1
   --twopassMode Basic"

   **⚪Input - sample.unmapped.fastq.gz `reads`, star (folder with star index) `index`. 🟢Output - sample.Aligned.sortedByCoord.out.bam `bam_sorted`, sample.Aligned.toTranscriptome.out.bam `bam_transcript`, sample.unmapped_1.fastq.gz `fastq`, sample.Log.final.out `log_final`, sample.Log.out `log_out`, sample.Log.progress.out `log_progress`**

   i. _SAMTOOLS_INDEX_GENOME_ : Indexes genome BAM file [Manual.](http://www.htslib.org/doc/samtools-index.html)

   **⚪Input - sample.Aligned.sortedByCoord.out.bam `bam`, 🟢Output - sample.Aligned.sortedByCoord.out.bam.bai `bai`.**

   ii. _SAMTOOLS_SORT_TRANSCRIPT_ : Sorts transcript bam file. [Manual.](http://www.htslib.org/doc/samtools-sort.html)

   **⚪Input - sample.Aligned.toTranscriptome.out.bam `bam`, 🟢Output - sample.Aligned.toTranscriptome.out.bam.bai `bai`.**

   iii. _SAMTOOLS_INDEX_TRANSCRIPT_ : Indexes transcript bam file. [Manual.](http://www.htslib.org/doc/samtools-index.html)

   **⚪Input - sample.Aligned.sortedByCoord.out.bam `bam`, 🟢Output - sample.Aligned.sortedByCoord.out.bam.bai `bai`.**

5. <u>FILTER_TRANSCRIPTS</u>
6. <u>BAM_SORT_STATS_SAMTOOLS_TRANSCRIPT</u>
7. <u>GENOME_DEDUP</u>
8. <u>TRANSCRIPT_DEDUP</u>
9. <u>CALC_GENOME_CROSSLINKS</u>
10. <u>CALC_TRANSCRIPT_CROSSLINKS</u>
11. _CLIPPY_
12. <u>PARACLU_ANALYSE</u>
13. <u>ICOUNT_ANALYSE</u>
14. _PEKA_
15. <u>CLIPSEQ_CLIPQC</u>
16. _MULTIQC_

## Common Issues

**Analysis of publicly available data**

Analysis of publicly available data, especially older data can present issues as sometimes the barcode structure or structure of the uploaded fastq is not obvious. Things to check are:

- Look at the first few reads in the fastq - is the random barcode (UMI) in the header? This tells us that the 5' barcode and 3' adapter have been removed from the read already, make sure to check what format the random barcode (UMI) is provided in, eg. "rbc:" and change the `umi_separator` parameter if you need to.
- If the random barcode is not in the header then you can use the `run_move_umi_to_header = true` parameter. The structure of the barcodes should be provided in the manuscript, but this is not always the case.
- Be careful for data in SRA, where the read header is entirely removed - this might have removed the information of random barcode without the manuscript authors noticing.

## FAQ

**What is the difference between this pipeline and nf-core/clipseq?**

This pipeline is based on nf-core/clipseq but is updated in a few ways. For one, the Nextflow code has been updated from DSL1 to the modular DSL2. This pipeline also includes additional tools, namely [PEKA](https://github.com/ulelab/peka) for positional motif analysis and [Clippy](https://github.com/ulelab/clippy) for peak calling.

**I would like to contribute to this pipeline, or make a request for a new feature. Where can I do this?**

## Authors

This DSL2 Nextflow pipeline is maintained by Goodwright. It was updated from the DSL1 [nf-core/clipseq](https://nf-co.re/clipseq) pipeline in collaboration with the original authors and Prof. Jernej Ule.

## References
