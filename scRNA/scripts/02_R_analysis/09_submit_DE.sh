#!/bin/sh
#SBATCH --job-name=DE                    # Job name
#SBATCH --mem=500G                       # Memory allocation
#SBATCH --output="logs/%x.%j.stdout"     # Standard output log file
#SBATCH --error="logs/%x.%j.stderr"      # Standard error log file
#SBATCH --time=7-00:00:00                # Time limit (7 days)

# Source environment variables (MY_DIR var)
source ../../refs/.env

# Load Singularity module
module load singularity

# Define Singularity image
SIMG_FILE_NAME="rstudio-4.3.0-4-with_modules.sif"
SIMG_IMAGE="/packages/containers/RStudio/${SIMG_FILE_NAME}"

# Ensure the same R library paths as in RStudio
# NOTE: MY_DIR is an environment variable I sourced at the start
export RPRIHOME="${MY_DIR}/R/${SIMG_FILE_NAME}/"
export R_LIBS_USER="${RPRIHOME}libs"
export HOME="${MY_DIR}"  

# Set up Singularity bind paths to match RStudio environment
BIND_PATHS="/tgen_labs,/opt,${RPRIHOME}:${RPRIHOME},${HOME}/.Rprofile:${HOME}/.Rprofile"

# Run the script inside Singularity with the correct environment
singularity exec -H "${MY_DIR}" --bind "${BIND_PATHS}" "${SIMG_IMAGE}" Rscript "08b_DE.R"

echo "DE job completed"
