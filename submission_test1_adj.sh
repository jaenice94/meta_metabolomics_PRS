#!/bin/bash
#SBATCH --job-name=PRScs
#SBATCH --output=logs/PRScs.%j.log
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G


cd /dss/dsshome1/0A/di54jot/meta_metabolomics/prscs
mkdir -p logs

################################################
# Define sample sheet, config file and snakemake files to use 
################################################

SAMPLESHEET="csv/gwas_list_test1.csv"
SNAKEFILE="${SNAKEFILE:-workflow/snakeflow_assoc.smk}" #keep
CONFIGFILE="${CONFIGFILE:-config/config_assoc.yaml}" #fill config with your paths


if [[ ! -f "$SAMPLESHEET" ]]; then
  echo "ERROR: Samplesheet does not exist: $SAMPLESHEET"
  exit 1
fi


#edit to your cluster requirements

# Load conda
module load miniforge3
source "$(conda info --base)/etc/profile.d/conda.sh"

conda activate snakemake7

#test dry run: snakemake -np -s "$SNAKEFILE" --configfile "$CONFIGFILE" --config "samplesheet=$SAMPLESHEET" 

snakemake \
  -s "$SNAKEFILE" \
  --configfile "$CONFIGFILE" \
  --config "samplesheet=$SAMPLESHEET" \
  --use-conda \
  --conda-frontend conda \
  --latency-wait 60 \
  --cores 1 \
  --keep-incomplete \
  --printshellcmds \
  --rerun-incomplete
