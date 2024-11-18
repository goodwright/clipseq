#!/bin/bash
#SBATCH --job-name=API_upload
#SBATCH --time=12:00:00  # Set the maximum runtime (e.g., 12 hour)
#SBATCH --ntasks=1  # Number of tasks (usually 1 for a single script)
#SBATCH --cpus-per-task=1  # Number of CPU cores per task
#SBATCH --mem=6G  # Memory per node (e.g., 4GB)
#SBATCH --nodes=1
#SBATCH --output=/nemo/lab/ulej/home/users/baih/PM24162/goodwright_clipseq/clipseq/tests/ARTR-seq/logs/API_upload_20240924_%j.out


# Run the Python script
ml Anaconda3/2023.09-0
conda activate flow
python3 /nemo/lab/ulej/home/users/baih/PM24162/goodwright_clipseq/clipseq/tests/ARTR-seq/script/API_upload.py