# PRS-Metabolite Association Pipeline

A snakemake wrapped pipeline to generate polygenic risk scores (PRS) using GWAS and test associations between PRS and metabolite traits. PRS are computed using PRS-CS and scored with PLINK1.9. Association models fit lipid ~ PRS + covariates (PRS-model) and lipid ~ covariates (baseline-model), reporting effect sizes (BETA), p-values and incremental R2 (PRS-model - baseline model). 

## Software requirements:
- snakemake 
- prs-cs
- plink1.9

All necessary software can be obtained as follows: 

Clone this git and the original PRScs.git
```bash
git clone https://github.com/jaenice94/prscs.git #modify!! 
cd prscs
git clone https://github.com/getian107/PRScs.git
```

## Data requirements:
- Formatted GWAS study (provide specifics) not containing the genotype to be tested (leave-one-out if necessary).
- LD reference (use ldblk_1kg_eur) - can be derived from XX. 
- Genotype data (you will need a .bim containing all relevant SNPs, .fam containing relevant samples, a list of .dosage files - specified in dosage_list.txt). The genotype data should be QC'd. Variants with MAF etc removed. 
- Normalised metabolite data (IID Peak1 Peak2 Peak3 ... ) measured in same individuals as the genotype. Please remove non-annotated lipids, fasting lipids and outliers. 
- Metadata containing information on covariates of same individuals as genotype and metabolites 

Before running the pipeline, edit the dosage_list.txt, config.yaml and samplesheet with your cohort specific information. 
- dosage_list.txt contains a list of paths to genotype dosage files *.dosage.gz from processing with plink19 
- gwas_list.csv please provide here the dir_path of the formatted GWAS and the N of the leave-one-out GWAS if applicable. 
- config.yaml please provide your directory paths for the necessary inputs

## Running the pipeline

Step 1 - set up environment 
```bash
conda env create -f workflow/envs/snakemake7.yaml
conda activate snakemake7
```
Step 2 - test out the snakemake.smk and your inputs in a dry run 

```bash
cd prscs 

SAMPLESHEET="csv/gwas_list_test1.csv"
SNAKEFILE="${SNAKEFILE:-workflow/snakeflow_assoc.smk}" #keep
CONFIGFILE="${CONFIGFILE:-config/config_assoc.yaml}" #fill config with your paths

snakemake -np -s "$SNAKEFILE" --configfile "$CONFIGFILE" --config "samplesheet=$SAMPLESHEET" 

Step 3A - submit pipeline as script (modify the #SBATCH header to fot your cluster's requirements on job submissions)

sbatch submit_assoc.sh 
```

or 

Step 3B - run from within tmux session (not tested; only from login node) - this is for slurm on LRZ; modify to fit your HPC requirements
```bash
cd prscs 

SAMPLESHEET="csv/gwas_list_test1.csv"
SNAKEFILE="${SNAKEFILE:-workflow/snakeflow_assoc_tmux.smk}" #snakefile specifying resources / or always specify and just have no effect in submit_assoc.sh
CONFIGFILE="${CONFIGFILE:-config/config_assoc.yaml}" #fill config with your paths

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
  --rerun-incomplete \
  --cluster "sbatch --clusters=serial --partition=serial_std --time={resources.runtime} --mem={resources.mem_mb}M --output=logs/{rule}.%j.log"
```
## Structure of repository
```
├──workflow/snakemake_assoc.smk #snakemake script
│  └──envs #enviroment.yamls required
│  └──scripts/prs_metabolite_associations.R #contains necessary R-scripts
├──PRScs/ #PRScs scripts cloned from https://github.com/getian107/PRScs
├──config/config_assoc.yaml #config file  - please enter your specifics
├──csv/gwas_list.csv #sheet for GWAS to run - please update your specifics
├──dosage_list.txt #paths to per-chromosome genotype dosage files - please update your specifics
├──results #default output directory for association results
└──submit_assoc.sh #SLURM submission script - please update to match your HPC
```
