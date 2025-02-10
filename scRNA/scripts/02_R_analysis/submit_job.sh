#!/bin/sh
#SBATCH --job-name=DE                    # Job name
#SBATCH --mem=500G                       # Memory allocation
#SBATCH --output=logs/%x.%j.stdout       # Standard output log file
#SBATCH --error=logs/%x.%j.stderr        # Standard error log file
#SBATCH --time=2-00:00:00                # Time limit (DD-HH:MM:SS)

# Load Singularity module
module load singularity

# Define Singularity image
SIMG_FILE_NAME=rstudio-4.3.0-4-with_modules.sif
if [ -d /packages/containers/RStudio/ ]; then
  SIMG_IMAGE=/packages/containers/RStudio/$SIMG_FILE_NAME
elif [ -d /opt/RStudio/ ]; then
  SIMG_IMAGE=/opt/RStudio/$SIMG_FILE_NAME
else
  echo "Singularity image not found!" && exit 1
fi

# Define script path
SCRIPT_PATH=07b_DE.R

# Define optional bind paths
BIND_PATHS="/tgen_labs,/opt"

# Run the script inside the Singularity container
singularity exec --bind $BIND_PATHS $SIMG_IMAGE Rscript $SCRIPT_PATH

echo "Script execution completed."
