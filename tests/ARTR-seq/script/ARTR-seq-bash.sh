#!/bin/sh
#SBATCH --job-name="ARTRseq_HeLa_results"
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --output=/nemo/lab/ulej/home/users/baih/PM24162/goodwright_clipseq/v1.3/clipseq/tests/ARTR-seq_20240925/logs/1st_run_20240925-%A.out

ml purge
ml Anaconda3/2023.09-0
ml Nextflow/23.04.4
ml Singularity/3.11.3
# ml Graphviz/2.38.0-foss-2016b

# the pipeline also has a "libraryDir = '/flask/apps/containers/all-singularity-images'" in the config file, which searches for the singularity images on nemo
export NXF_SINGULARITY_CACHEDIR=/camp/lab/ulej/home/shared/singularity
cd /nemo/lab/ulej/home/users/baih/PM24162/goodwright_clipseq/v1.3/clipseq

# Input Parameters
input_dir="/nemo/stp/sequencing/inputs/instruments/fastq/20240807_LH00442_0034_A22MYGKLT3/fastq/PM24162"
outdir="/nemo/lab/ulej/home/users/baih/PM24162/goodwright_clipseq/v1.3/clipseq/tests/ARTR-seq_20240925/output"
sample_sheet=/nemo/lab/ulej/home/users/baih/PM24162/goodwright_clipseq/v1.3/clipseq/tests/ARTR-seq_20240925/data/sample_sheet.csv

genome_fasta="/nemo/lab/ulej/home/users/baih/Reference/Gencode_reference_files/Human/GRCh38.v44/GRCh38.primary_assembly.genome.fa"
smRNA_fasta="/nemo/lab/ulej/home/users/baih/Reference/Gencode_reference_files/Human/GRCh38.v44/flow_bio_additional_ref/Homo_sapiens.GRCh38.smrna.fasta"
genome_gtf="/nemo/lab/ulej/home/users/baih/Reference/Gencode_reference_files/Human/GRCh38.v44/gencode.v44.primary_assembly.annotation.gtf"

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
umi_header_format="NNNNN"
umi_separator="_"

nextflow run main.nf -resume \
-profile crick \
-N pelusobai@gmail.com \
--outdir $outdir \
--samplesheet $sample_sheet \
--fasta $genome_fasta \
--smrna_fasta $smRNA_fasta \
--gtf $genome_gtf \
--move_umi_to_header true \
--umi_header_format "${umi_header_format}" \
--trimgalore_params "${trimgalore_params}" \
--bowtie_params "${bowtie_params}" \
--star_params "${star_params}" \
--umi_separator "${umi_separator}" 