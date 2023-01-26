# CLIP Analysis

Full documentation [here](Documentation.md), salient usage details summarised below.

Cross-linking and immunoprecipitation followed by sequencing (CLIP) has allowed high resolution studies of RNA binding protein (RBP)-RNA interactions at transcriptomic scale.

This pipeline enables analysis of various forms of single-end CLIP data including variants of iCLIP (eg. irCLIP, iCLIP2, iiCLIP) and eCLIP (note that we have achieved comparable results using this pipeline to study reformatted paired-end eCLIP also). The pipeline currently doesn't support mutation calling and therefore might not be suitable for PAR-CLIP analysis, but we plan to include this at a future date.

## Quick-start

To test the pipeline, use the associated config file and run it with the
profiles `test` and the container engine you wish to use eg. `docker`. For example:

```bash
nextflow run main.nf -profile test,docker
```

Full dataset testing of 9 iCLIP samples can also be run using profile `test_full`.
A test can also be run that skips all preparing of annotations/indexes using profile `test_no_prep_genome`.

## Input

If you require all reference files (eg. genomic indexes, filtered and segmented gtf...) to be generated the minimal input is:

- `samplesheet` : csv file containing 4 columns: group,replicate,fastq_1,fastq_2. group is the sample name, replicate is currently unused by the pipeline so filling with '1' is acceptable, fastq_1 is your demultiplexed sample fastq, paired end is currently not supported so please do not add a fastq 2 .eg './tests/data/samplesheets/small-single-sample-se.csv'

| group      | replicate |  fastq_1                                                                | fastq_2 |
| ---------- | --------- | ----------------------------------------------------------------------  | ------- |
| TDP43_1    | 1         | s3://nf-core-awsmegatests/clipseq/input_data/fastq/ERR1530360.fastq.gz  |         |


- `fasta`       : genome fasta file .eg './tests/data/genome/homosapien-hg37-chr21.fa.gz'
- `smrna_fasta` : fasta file to be mapped to before the genome file, typically containing rRNA and tRNA sequences  .eg'./tests/data/genome/homosapiens_smallRNA.fa.gz'
- `gtf`         : annotation file for the genome fasta .eg'./tests/data/genome/gencode.v35.chr21.gtf.gz'

If you are providing all reference files then the following *additional* files must be provided (note these are all produced by the `prepare_clipseq` subworkflow:

- `fasta_fai`
- `chrom_sizes`
- `genome_index`
- `smrna_index`
- `smrna_fasta_fai`
- `smrna_chrom_sizes`
- `longest_transcript`
- `filtered_gtf`
- `seg_gtf`
- `seg_filt_gtf`
- `seg_resolved_gtf`
- `seg_resolved_gtf_genic`
- `regions_gtf`
- `regions_filt_gtf`
- `regions_resolved_gtf`
- `regions_resolved_gtf_genic`

## Output

When the full pipeline is run, output is organised into 6 folders:
- `00_genome` contains all reference files produced when the prepare_clipseq subworkflow is run.
- `01_prealign` contains pre-trimmed FastQC reports, trimmed read files and post-trimming FastQC reports.
- `02_alignment` contains two folders, one for the pre-mapping "smrna" and one for the genomic mapping "target", each contain alignment files and the "target" folder also contains useful samtools assessment of the bam file.
- `03_filt_dedup` contains de-duplicated genome mapped bams along with statistics and also two transcriptome mapped de-duplicated bams. Here those marked "filt" have been filtered to only contain alignments to the longest transcript for each gene.
- `04_crosslinks` contains genomic crosslink bed, bedgraph and normalised bedgraph (crosslinks divided by total number of crosslinks in the sample and multiplied by a million resulting in a crosslinks per million (CPM) value); also the same files for transcriptome mapping filtered by longest transcript.
- `05_peak_calling` contains Clippy, iCount and Paraclu peaks. The iCount folder also contains gene, subtype and type level summaries of crosslink information and metagene plots around transcript landmarks of interest in the rnamaps folder. Also included is PEKA output.
- `06_reports` contains various CLIP-specific QC metrics in tabular format in the clipqc folder. These are plotted, alongside other QC metrics in the html provided in the multiqc folder.

## Authors and contact
This DSL2 Nextflow pipeline is maintained by Goodwright. It was updated from the DSL1 nf-core/clipseq pipeline in collaboration with the original authors and Prof. Jernej Ule. 
To raise any issues or comments with the pipeline you can (in order of preference):
- Raise an issue in this repository
- Write to us in our [Slack](https://join.slack.com/t/imapsgroup/shared_invite/zt-r24y3591-Xbhnym2t38u_urU~I0K0lQ) 
- Email charlotte.capitanchik@goodwright.com
