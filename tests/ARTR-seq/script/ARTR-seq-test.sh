#!/bin/sh
#SBATCH --job-name="ARTRseq_test"
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --output=/nemo/lab/ulej/home/users/baih/PM24162/goodwright_clipseq/clipseq/tests/ARTR-seq/logs/ARTR_seq_test_20241024-%A.out

ml purge
ml Anaconda3/2023.09-0
ml Nextflow/23.04.4
ml Singularity/3.11.3
# ml Graphviz/2.38.0-foss-2016b

# the pipeline also has a "libraryDir = '/flask/apps/containers/all-singularity-images'" in the config file, which searches for the singularity images on nemo
export NXF_SINGULARITY_CACHEDIR=/camp/lab/ulej/home/shared/singularity
cd /nemo/lab/ulej/home/users/baih/PM24162/goodwright_clipseq/clipseq

# Input Parameters
outdir="/nemo/lab/ulej/home/users/baih/PM24162/goodwright_clipseq/clipseq/tests/ARTR-seq/output"
sample_sheet=/nemo/lab/ulej/home/users/baih/PM24162/goodwright_clipseq/clipseq/tests/ARTR-seq/input_data/sample_sheet_for_testing.csv
pairwise_samplesheet=/nemo/lab/ulej/home/users/baih/PM24162/goodwright_clipseq/clipseq/tests/ARTR-seq/input_data/pairwise_peak_calling_for_testing.csv

# Reference Parameters
fasta="/nemo/lab/ulej/home/users/baih/Reference/Gencode_reference_files/Human/GRCh38.v44/GRCh38.primary_assembly.genome.fa"
smRNA_fasta="/nemo/lab/ulej/home/users/baih/Reference/Gencode_reference_files/Human/GRCh38.v44/flow_bio_additional_ref/Homo_sapiens.GRCh38.smrna.fasta"
genome_gtf="/nemo/lab/ulej/home/users/baih/Reference/Gencode_reference_files/Human/GRCh38.v44/gencode.v44.primary_assembly.annotation.gtf"

target_genome_index=/nemo/lab/ulej/home/users/baih/Reference/clip-seq_prepared_genomes/00_genome/index/star
smrna_genome_index=/nemo/lab/ulej/home/users/baih/Reference/clip-seq_prepared_genomes/00_genome/index/bowtie
fasta_fai=/nemo/lab/ulej/home/users/baih/Reference/clip-seq_prepared_genomes/00_genome/GRCh38.primary_assembly.genome.fa.fai
filtered_gtf=/nemo/lab/ulej/home/users/baih/Reference/clip-seq_prepared_genomes/00_genome/gencode_filtered.gtf
chrom_sizes=/nemo/lab/ulej/home/users/baih/Reference/clip-seq_prepared_genomes/00_genome/GRCh38.primary_assembly.genome.fa.sizes
smrna_fasta_fai=/nemo/lab/ulej/home/users/baih/Reference/clip-seq_prepared_genomes/00_genome/Homo_sapiens.GRCh38.smrna.fasta.fai
smrna_chrom_sizes=/nemo/lab/ulej/home/users/baih/Reference/clip-seq_prepared_genomes/00_genome/Homo_sapiens.GRCh38.smrna.fasta.sizes
longest_transcript=/nemo/lab/ulej/home/users/baih/Reference/clip-seq_prepared_genomes/00_genome/longest_transcript.txt
longest_transcript_fai=/nemo/lab/ulej/home/users/baih/Reference/clip-seq_prepared_genomes/00_genome/longest_transcript.fai
longest_transcript_gtf=/nemo/lab/ulej/home/users/baih/Reference/clip-seq_prepared_genomes/00_genome/longest_transcript.gtf
seg_gtf=/nemo/lab/ulej/home/users/baih/Reference/clip-seq_prepared_genomes/00_genome/gencode_seg.gtf
seg_filt_gtf=/nemo/lab/ulej/home/users/baih/Reference/clip-seq_prepared_genomes/00_genome/gencode_filtered_seg.gtf
seg_resolved_gtf=/nemo/lab/ulej/home/users/baih/Reference/clip-seq_prepared_genomes/00_genome/gencode_filtered_seg_genicOtherfalse.resolved.gtf
seg_resolved_gtf_genic=/nemo/lab/ulej/home/users/baih/Reference/clip-seq_prepared_genomes/00_genome/gencode_filtered_seg_genicOthertrue.resolved.gtf
regions_gtf=/nemo/lab/ulej/home/users/baih/Reference/clip-seq_prepared_genomes/00_genome/gencode_regions.gtf
regions_filt_gtf=/nemo/lab/ulej/home/users/baih/Reference/clip-seq_prepared_genomes/00_genome/gencode_filtered_regions.gtf.gz
regions_resolved_gtf=/nemo/lab/ulej/home/users/baih/Reference/clip-seq_prepared_genomes/00_genome/gencode_filtered_regions_genicOtherfalse.resolved.gtf
regions_resolved_gtf_genic=/nemo/lab/ulej/home/users/baih/Reference/clip-seq_prepared_genomes/00_genome/gencode_filtered_regions_genicOthertrue.resolved.gtf

# Pipeline Parameters
# here we remove 3 nucleotides from the 5' end of read 1 and the 3' end of read 2, which are bias introduced by MMLV
# also remove 4 nucleotides from the 3' end of read 1 and the 5' end of read 2, which are bias introduced by random RT priming
trimgalore_params="--fastqc --length 10 -q 20 --clip_r1 3 --three_prime_clip_r1 4"
# Here it's the bowtie2 alignment parameters, basically we allow 2 mismatches and 100 alignments, and we only keep the best alignment
# Note that the --best and --strata do not work in pared-end mode, but here we are only interested in the unaligned reads so it's fine
# Also note that ARTR-seq used bowtie2 -n mode with seed legnth 15 ()
bowtie_params="-v 2 -m 100 --norc --best --strata"
# Here we use the STAR parameters in CLIP-seq, which allow multimappers and basically a classic encode RNA-seq pipeline.
# Note the alignEndsType is EndToEnd, which prevents soft-clipping, since we hard-clipped the reads in the previous step
star_params="--outFilterMultimapNmax 100 --outFilterMultimapScoreRange 1 --outSAMattributes All --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterType BySJout --alignIntronMin 20 --alignIntronMax 1000000 --outFilterScoreMin 10 --alignEndsType EndToEnd --twopassMode Basic"
macs3_params="--keep-dup all --nomodel --extsize 30 --bdg -g hs"
umi_header_format="NNNNN"
umi_separator="_"

nextflow run main.nf -resume \
-profile crick \
-N pelusobai@gmail.com \
--whole_read_analysis true \
--outdir $outdir \
--samplesheet $sample_sheet \
--pairwise_samplesheet $pairwise_samplesheet \
--fasta $fasta \
--smrna_fasta $smRNA_fasta \
--gtf $genome_gtf \
--target_genome_index $target_genome_index \
--smrna_genome_index $smrna_genome_index \
--fasta_fai $fasta_fai \
--filtered_gtf $filtered_gtf \
--chrom_sizes $chrom_sizes \
--smrna_fasta_fai $smrna_fasta_fai \
--smrna_chrom_sizes $smrna_chrom_sizes \
--longest_transcript $longest_transcript \
--longest_transcript_fai $longest_transcript_fai \
--longest_transcript_gtf $longest_transcript_gtf \
--seg_gtf $seg_gtf \
--seg_filt_gtf $seg_filt_gtf \
--seg_resolved_gtf $seg_resolved_gtf \
--seg_resolved_gtf_genic $seg_resolved_gtf_genic \
--regions_gtf $regions_gtf \
--regions_filt_gtf $regions_filt_gtf \
--regions_resolved_gtf $regions_resolved_gtf \
--regions_resolved_gtf_genic $regions_resolved_gtf_genic \
--move_umi_to_header true \
--umi_header_format "${umi_header_format}" \
--trimgalore_params "${trimgalore_params}" \
--bowtie_params "${bowtie_params}" \
--star_params "${star_params}" \
--macs3_params "${macs3_params}" \
--umi_separator "${umi_separator}"