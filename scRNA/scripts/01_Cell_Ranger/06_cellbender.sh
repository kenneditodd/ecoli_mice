#!/bin/sh
#SBATCH --job-name=cellbender       # Name of the job
#SBATCH --mem=50G                   # Amount of memory allocated for the job
#SBATCH --output logs/%x.%j.stdout  # File for standard output
#SBATCH --error=logs/%x.%j.stderr   # File for standard error output
#SBATCH --partition=gpu             # Specifies the partition (queue)
#SBATCH --gres=gpu:1                # Specify GPU
#SBATCH --cpus-per-task=16          # Number of CPU cores per task
#SBATCH --time=24:00:00             # Maximum time the job is allowed to run, HH:MM:SS

# source settings
source $HOME/.bash_profile
source ../../refs/.env

# activate conda environment
conda activate cellbender

# cellbender version
cellbender --version

# print sample passed from sample_loop.sh script
SAMPLE=$1
echo "sample: $SAMPLE"

# in order for ckpt.tar.gz in cellbender to save correctly, go to where data is located
cd ${PROJECT_DIR}/counts/${SAMPLE}/outs

# run cellbender remove-background
cellbender remove-background \
      --input raw_feature_bc_matrix.h5 \
      --output ${SAMPLE}_cellbender.h5 \
      --cuda

# Key
# --input       :un-filtered data file, .h5 or matrix directory from 10x
# --output      :output file location containing .h5 extension
# --cuda        :for use with GPU
