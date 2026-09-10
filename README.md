# PRS-Metabolite Association Pipeline

A snakemake workflow to generate polygenic risk scores (PRS) from GWAS summary statistics using PRS-CS, score PRS in a target cohort using PLINK2 and test associations between the PRS and metabolic traits. 

### Pipeline workflow: 
<img width="1400" height="788" alt="Image" src="https://github.com/user-attachments/assets/3a8f5eb7-6bbc-45f1-bbbe-7400293d8b1c" />
Pipeline overview figure generated with ChatGPT.

## 1. Set-up:

### Software requirements:

The workflow uses:
```
- Snakemake 7 
- PRS-CS
- PLINK2
- Python 3.10
- R
- conda
```
Please ensure that conda is available on your system. The software required by individual workflow steps is specified in the YAML enviroment under workflow/envs/ and is installed automatically by Snakemake during the workflow. 


### Clone Github Repositories:

clone this repository and inside it the original PRS-CS repository:
```bash
git clone https://github.com/jaenice94/meta_metabolomics_PRS.git
cd meta_metabolomics_PRS
git clone https://github.com/getian107/PRScs.git
```
### Create Snakemake environment: 
```bash
conda env create -f workflow/envs/snakemake7.yaml
conda activate snakemake7
```

## 2. Data requirements:

### GWAS summary statistics:
each GWAS must contain SNP identifier, effect allele, other allele, effect estimate (BETA or OR) and either standard error (SE) or p-value (P). If available BETA + SE are preferred. The original column names can be specified in csv/gwas_list.csv, the workflow then extracts and renames columns for PRS-CS.

Note: avoid sample overlap between GWAS and target cohort. For cohorts contributing to the original GWAS/meta-analysis, use leave-one-cohort-out GWAS summary statistics and specify GWAS sample size in csv/gwas_list.csv. 

### LD reference: 
Download the appropriate PRS-CS LD reference (see https://github.com/getian107/PRScs.git for download links) - use ldblk_1kg_eur for PGC meta-metabolomics analysis unless otherwise specified.

Note: The LD reference is used internally by PRS-CS and to restrict and harmonise cohort variants to the HapMap3 SNP representation used by PRS-CS. Matching is done by rsID, chromosome and allele and should therefore work for both genotypes in GRCh37 and GRCh38. 

### Genotype data: 
The pipeline currently supports PLINK1 (.bed,.bim,.fam) and PLINK2 (.pgen, .pvar, .psam) files that are chromosome-split or genome-wide. For PLINK2 both best-guess and dosage information genotypes are supported. Target cohort genotype data must have undergone appropriate QC post-imputation. The sample IID must match to sample IID in metabolite and metadata files.

### Normalised metabolite data: 
For meta-metabolomics analayis, use the normalisd metabolite output.

Note: Must contain IID and allows also a FID column. All other columns are treated as metabolite traits and will enter the association analysis. Ensure missing values are encoded as NA. 

### Metadata: 
Metadata file contains IID and all covariates to be used in the PRS-association model, specified in config/config_cohort.yaml.

Note: Sample IID must correspond to IID in genotype and metabolite files. Ensure missing values are encoded as NA. 
  
## 3. Configure pipeline files: 

### GWAS samplesheet:

csv/gwas_list.csv is a ;-separated file that contains per row one GWAS phenotype for which a PRS will be calculated. See csv/gwas_list.csv for example and column input information. 

### Cohort configuration:

config/config_cohort.yaml contains cohort and analysis-specific configurations. Please update file to include your cohort-specific inputs. 

## 4. Test configurations with dry-run: 

From within repository directory:
```{bash}

cd prscs

#create log directory
mkdir -p logs

conda activate snakemake7

################################################
# Define sample sheet, config file and snakemake files to use 
################################################

SAMPLESHEET="csv/gwas_list.csv"
CONFIGFILE="config/config_cohort.yaml"

SNAKEFILE="workflow/snakemake_harmonised.smk"

snakemake -np -s "$SNAKEFILE" --configfile "$CONFIGFILE" --config "samplesheet=$SAMPLESHEET" 
```

## 5. Run the pipeline. 

### A) with job submission (here an example script for SLURM): 
```{bash}
#!/bin/bash
#SBATCH --job-name=PRScs
#SBATCH --output=logs/PRScs.%j.log
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --time=16:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

################################################
# Define sample sheet, config file and snakemake files to use 
################################################

SAMPLESHEET="csv/gwas_list.csv" #updated with your paths
CONFIGFILE="config/config_cohort.yaml" #updated with your paths

SNAKEFILE="workflow/snakemake_harmonised.smk" #keep

if [[ ! -f "$SAMPLESHEET" ]]; then
  echo "ERROR: Samplesheet does not exist: $SAMPLESHEET"
  exit 1
fi

echo "$SAMPLESHEET"
echo "$SNAKEFILE"
echo "$CONFIGFILE"

module load miniforge3
conda activate snakemake7

snakemake \
  -s "$SNAKEFILE" \
  --configfile "$CONFIGFILE" \
  --config "samplesheet=$SAMPLESHEET" \
  --use-conda \
  --conda-frontend conda \
  --latency-wait 60 \
  --cores 8 \
  --keep-incomplete \
  --printshellcmds \
  --rerun-incomplete
```
Submit from within repository directory with: 
```{bash}
sbatch submit_snakemake.sh
```

### B) Run within a tmux session (not extensively tested): 


### Unlocking an interrupted Snakemake run: 

If the previous snakemake run was interrupted, Snakemake will report that the working directory is locked. Unlock it with: 
```{bash}
SAMPLESHEET="csv/gwas_list.csv" 
CONFIGFILE="config/config_cohort.yaml" 
SNAKEFILE="workflow/snakemake_harmonised.smk" 

snakemake -s "$SNAKEFILE" --configfile "$CONFIGFILE" --config "samplesheet=$SAMPLESHEET" --unlock --cores 1
```
## Repository structure
```
├──workflow/snakemake_harmonised.smk #snakemake pipeline
│  └──envs #enviroment.yamls required
│        └──prscs.yaml
│        └──r_base.yaml
│        └──plink.yaml
│        └──snakemake7.yaml
│  └──scripts #contains necessary R-scripts
│        └──format_GWAS.R
│        └──harmonise_bim_hm3.R
│        └──sum_scores_p12.R
│        └──prs_metabolite_associations_p12.R
├──PRScs/ #PRScs scripts cloned from https://github.com/getian107/PRScs
│        └──PRScs.py
├──config #config files
│        └──config_cohort.yaml #update to your cohort-specifics 
├──csv #samplesheet for GWAS to process
│        └──gwas_list.csv #update to your specifics 
└──submit_snakemake.sh #example of submission script for SLURM
└──FADS_variants_GRCh37.txt #FADS variants to exclude
└──exclude_variants.sh #example to derive variant to exclude
```
