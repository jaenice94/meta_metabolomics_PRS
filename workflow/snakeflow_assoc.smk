############################################################
#read in configs from CSV
############################################################
import pandas as pd
from pathlib import Path

SAMPLESHEET = config.get("samplesheet") 				#which is provided via snakemake --config samplesheet=path/to/file.csv
SEP = config.get("base", {}).get("sep", ";")

#throw error if no samplesheet
if not SAMPLESHEET:
	raise ValueError("No samplesheet provided. Pass --config samplesheet=PATH")

#read samplesheet
df = pd.read_csv(SAMPLESHEET, sep=SEP, dtype=str, comment="#").fillna("")

#specify required columns
required = ["phenotype","sst_file","n_gwas"]
missing = [c for c in required if c not in df.columns]
if missing:
	raise ValueError(f"Samplesheet missing columns: {missing}")

############################################################
#read in general settings from config file
############################################################

#general settings
general = config.get("general", {})
ld_ref = general.get("ld_ref_dir")
bim = general.get("bim")
genotype = general.get("genotype")
list_file = general.get("list_file")
chromosomes = general.get("chromosomes")
fam=general.get("fam")
output_dir = general.get("output_dir", "results")
association = config.get("association_analysis", {})

############################################################
#Provide definitions
############################################################

#get run parameters
RUNS = df.to_dict("records")

RUN_BY_PHENO = {r["phenotype"]: r for r in RUNS}
def get_run(wc):
    return RUN_BY_PHENO[wc.phenotype]

############################################################
#Define the final outputs
############################################################

rule all:
	input:
		expand(
			f"{output_dir}/{genotype}_{{phenotype}}_phi_auto/{genotype}_{{phenotype}}_phi_auto.profile",
			phenotype=df["phenotype"].tolist()
		),
		f"results/PRS_metabolite_associations_{genotype}.txt",


############################################################
#Run PRSCS - here one run per phenotype/prscs settings
############################################################

rule run_prscs:
	input:
		ld_ref=ld_ref,
		bim=f"{bim}.bim", 
		gwas=lambda wc: get_run(wc)["sst_file"]
	output:
		log=f"{output_dir}/{genotype}_{{phenotype}}_phi_auto/{genotype}_{{phenotype}}_phi_auto_pst_eff_a1_b0.5_phiauto_chr{{chr}}.log",
		weights=f"{output_dir}/{genotype}_{{phenotype}}_phi_auto/{genotype}_{{phenotype}}_phi_auto_pst_eff_a1_b0.5_phiauto_chr{{chr}}.txt",
	params:
		bim=bim,
		prscs=f"PRScs/PRScs.py",
		n_gwas=lambda wc: get_run(wc)["n_gwas"],
		outdir=lambda wc: f"{output_dir}/{genotype}_{wc.phenotype}_phi_auto",
		prefix=lambda wc: f"{output_dir}/{genotype}_{wc.phenotype}_phi_auto/{genotype}_{wc.phenotype}_phi_auto"
	conda:
		"envs/prscs.yaml"
	shell:
		r"""
		
		mkdir -p {params.outdir}

		echo "Running PRSCS with {input.gwas} for chr{wildcards.chr}"

		python3 {params.prscs} \
  			--ref_dir={input.ld_ref} \
  			--bim_prefix={params.bim} \
  			--sst_file={input.gwas} \
  			--n_gwas={params.n_gwas} \
  			--seed=42 \
  			--chrom={wildcards.chr} \
  			--out_dir={params.prefix} \
  			> {output.log} 2>&1
		"""

#done - test! 

rule combine_weights:
	input:
		weights=lambda wc: expand(
			f"{output_dir}/{genotype}_{wc.phenotype}_phi_auto/{genotype}_{wc.phenotype}_phi_auto_pst_eff_a1_b0.5_phiauto_chr{{chr}}.txt",
			chr=chromosomes
			)
	output:
		weights=f"{output_dir}/{genotype}_{{phenotype}}_phi_auto/{genotype}_{{phenotype}}_phi_auto_scoring.txt"
	shell:
		r"""
		cat {input.weights} > {output.weights}
		"""

############################################################
#Calculate PRSCS for genotype
############################################################

rule score_PRS:
	input:
		dosage=list_file,
		fam=fam,
		weights=lambda wc: f"{output_dir}/{genotype}_{wc.phenotype}_phi_auto/{genotype}_{wc.phenotype}_phi_auto_scoring.txt",
	output:
		profile=f"{output_dir}/{genotype}_{{phenotype}}_phi_auto/{genotype}_{{phenotype}}_phi_auto.profile"
	params:
		profile=lambda wc: f"{output_dir}/{genotype}_{wc.phenotype}_phi_auto/{genotype}_{wc.phenotype}_phi_auto",
		genotype=genotype,

	conda:
		"envs/plink.yaml"
	shell:
		r"""
		set -euo pipefail

		echo "genotype:     {params.genotype}"
		echo "weights:  {input.weights}"
		echo "phenotype:  {wildcards.phenotype}"

		plink \
  			--dosage {input.dosage} list format=3 Zin \
  			--fam {input.fam} \
  			--score {input.weights} 2 4 6 no-mean-imputation \
  			--out {params.profile}

		echo  "calculating PRS for {wildcards.phenotype} for {params.genotype} done"
		"""

############################################################
#Run association models
############################################################

rule association_models: 
	input:
		profiles=expand(
			f"{output_dir}/{genotype}_{{phenotype}}_phi_auto/{genotype}_{{phenotype}}_phi_auto.profile",
			phenotype=df["phenotype"].tolist()
			),
		metabolites=association["metabolites"],
		metadata=association["metadata"]
	output:
		results=f"results/PRS_metabolite_associations_{genotype}.txt"
	params:
		covar=association["covariates"],
		phenotypes=df["phenotype"].tolist(),
	conda:
		"envs/r_base.yaml"
	script:
		"scripts/prs_metabolite_associations.R"

#check rule association models + rule all 
#check if we can harmonise covariate names or if we can otherwise harmonise the script 
