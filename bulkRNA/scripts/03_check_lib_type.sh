#!/bin/sh
#SBATCH --job-name=check_lib_type        # Job name
#SBATCH --mem=10G                        # Memory allocation
#SBATCH --output=logs/%x.%j.stdout       # Standard output log file
#SBATCH --error=logs/%x.%j.stderr        # Standard error log file
#SBATCH --tasks=15                       # Number of tasks (threads)
#SBATCH --time=02:00:00                  # Time limit (HH:MM:SS)

# Source environment variables
source ../refs/.env

# Load environment and activate the 'salmon' environment
source $HOME/.bash_profile
conda activate salmon

# Print the Salmon version
salmon -v

# Change to the raw reads directory
cd $RAW_READS_DIR

# Run Salmon quantification with validation of mappings
salmon quant --libType A \
             --index $SALMON_INDEX \
             --mates1 "${RAW_READS_DIR}/young_9_ms_t1_S9_L001_R1_001.fastq.gz" \
             --mates2 "${RAW_READS_DIR}/young_9_ms_t1_S9_L001_R2_001.fastq.gz" \
             --output "${PROJECT_DIR}/refs/transcript_quant" \
             --threads 15 \
             --validateMappings

# Key Information:
# --libType A : Autodetect library type
# --index : Path to the Salmon index
# --mates1 : FASTQ file for Read 1
# --mates2 : FASTQ file for Read 2
# --threads : Number of threads to use
# --validateMappings : Perform mapping validation

# Results:
# salmon (selective-alignment-based) v1.10.1
# Salmon will autodetect the library type, such as ISR (inward stranded reverse).
# ISR corresponds to '-s2' argument for featureCounts (reversely stranded).
