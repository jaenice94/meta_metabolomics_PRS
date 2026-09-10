############################################################
#read in configs from CSV
############################################################
import pandas as pd
from pathlib import Path

SAMPLESHEET = config.get("samplesheet")         #which is provided via snakemake --config samplesheet=path/to/file.csv
SEP = config.get("base", {}).get("sep", ";")

#throw error if no samplesheet
if not SAMPLESHEET:
  raise ValueError("No samplesheet provided. Pass --config samplesheet=PATH")

#read samplesheet
df = pd.read_csv(SAMPLESHEET, sep=SEP, dtype=str, comment="#").fillna("")

#specify required columns
required = ["phenotype","sst_file","n_gwas", "snp_col", "effect_allele_col", "other_allele_col", "effect_col", "effect_type", "stat_col", "stat_type"]
missing = [c for c in required if c not in df.columns]
if missing:
  raise ValueError(f"Samplesheet missing columns: {missing}")

valid_effect_types = {"BETA", "OR"}

invalid_effect = df.loc[
    ~df["effect_type"].str.upper().isin(valid_effect_types),
    ["phenotype", "effect_type"]
]

if not invalid_effect.empty:
    raise ValueError(
        "GWAS effect_type must be BETA or OR:\n"
        + invalid_effect.to_string(index=False)
    )

valid_stat_types = {"SE", "P"}

invalid_stat = df.loc[
    ~df["stat_type"].str.upper().isin(valid_stat_types),
    ["phenotype", "stat_type"]
]

if not invalid_stat.empty:
    raise ValueError(
        "GWAS Stat_col must be SE or P:\n"
        + invalid_stat.to_string(index=False)
    )

############################################################
#read in general settings from config file
############################################################

#general settings
general = config.get("general", {})
ld_ref = general.get("ld_ref_dir")
cohort = general.get("cohort")
gwas_dir = general.get("gwas_dir")
chromosomes = general.get("chromosomes")
output_dir = general.get("output_dir", "results")
Path(output_dir).mkdir(parents=True, exist_ok=True) #create if it does not exist yet

genotype_prefix = general.get("genotype_prefix")

association = config.get("association_analysis", {})
covariates_name = association.get("covariate_name")

hapmap_snps = str(Path(ld_ref) / "snpinfo_1kg_hm3")

region_exclusion = config.get("region_exclusion", {})
exclude_region = region_exclusion.get("exclude", False)
exclude_locus = region_exclusion.get("locus")
exclude_variants = region_exclusion.get("variants_to_exclude")

if exclude_region:
  if not exclude_locus:
    raise ValueError(
      "You need to define locus name for region_exclusion.locus"
      )

  if not exclude_variants:
    raise ValueError(
      "You need to provide a list of variants to exclude under region_exclusion.locus"
      )

  exclusion_suffix = f"excl_{exclude_locus}"

else: exclusion_suffix = ""


#validate genotype layout input
genotype_layout = general.get("genotype_layout")

if genotype_layout not in {"chromosome_split", "genomewide"}:
  raise ValueError(
    "general.genotype_layout must be "
    "'chromosome_split' or 'genomewide'"
    )

#validate genotype format input
genotype_format = general.get("genotype_format")
if genotype_format not in {"plink1", "plink2"}:
    raise ValueError(
        "general.genotype_format must be 'plink1' or 'plink2'"
    )

#validate naming structure - chr split files need chr indicator
if genotype_layout == "chromosome_split" and "{chr}" not in genotype_prefix:
  raise ValueError(
    "For chromosome_split genotype data, "
    "general.genotype_prefix must contain '{chr}'"
    )

if genotype_layout == "genomewide" and "{chr}" in genotype_prefix:
  raise ValueError(
    "For genomewide genotype data, "
    "general.genotype_prefix must not contain '{chr}'"
    )

############################################################
#Provide definitions
############################################################

#get run parameters
RUNS = df.to_dict("records")

RUN_BY_PHENO = {r["phenotype"]: r for r in RUNS}
def get_run(wc):
    return RUN_BY_PHENO[wc.phenotype]


def get_genotype_prefix(chrom):
  """
  Get chromosome-specific genotype prefix. 
  """
  return genotype_prefix.format(chr=chrom)


def get_variant_files():
  """
  Variant file(s) needed to create the PRS-CS BIM.
  """

  extension = ".pvar" if genotype_format == "plink2" else ".bim"

  if genotype_layout == "chromosome_split":
    return[
      f"{get_genotype_prefix(chrom)}{extension}"
      for chrom in chromosomes
    ]

  elif genotype_layout == "genomewide":
    return [
      f"{genotype_prefix}{extension}"
    ]

  else:
    raise ValueError(
      f"Unsupported genotype_layout: {genotype_layout}"
    )

def get_genotype_files(wc):
  """
  Genotype files required for scoring one chromosome.
  """

  prefix = get_genotype_prefix(wc.chr)

  if genotype_format == "plink2":
    return[
      f"{prefix}.pgen",
      f"{prefix}.pvar",
      f"{prefix}.psam"
      ]
  else:
    return[
      f"{prefix}.bim",
      f"{prefix}.bed",
      f"{prefix}.fam"
      ]

def get_genomewide_genotype_files(wc):
  """
  Genotype files required for genome-wide scoring.
  """

  prefix = genotype_prefix

  if genotype_format == "plink2":
    return[
      f"{prefix}.pgen",
      f"{prefix}.pvar",
      f"{prefix}.psam"
      ]
  else:
    return[
      f"{prefix}.bim",
      f"{prefix}.bed",
      f"{prefix}.fam"
      ]

#define which scoring pattern to use 
def get_score_files(wc):
  """
  Return score file(s) depending on genotype layout.
  """

  if genotype_layout == "chromosome_split":

    return expand(
      rules.score_PRS_chr.output.sscore,
      phenotype=wc.phenotype,
      chr=chromosomes
      )

  elif genotype_layout == "genomewide":

    return expand(
      rules.score_PRS_genomewide.output.sscore,
      phenotype=wc.phenotype
      )

#equivalent for scoring files with regional exclusion:
def get_excl_score_files(wc):
  """
  Return region-excluded score file(s) depending on genotype layout.
  """

  if genotype_layout == "chromosome_split":

    return expand(
      rules.score_PRS_chr_excl.output.sscore,
      phenotype=wc.phenotype,
      chr=chromosomes
      )

  elif genotype_layout == "genomewide":

    return expand(
      rules.score_PRS_genomewide_excl.output.sscore,
      phenotype=wc.phenotype
      )

plink_input_flag = (
  "--pfile"
  if genotype_format == "plink2"
  else "--bfile")

#bim files
prscs_bim_prefix = (
  f"{output_dir}/{cohort}_prscs"
  )
prscs_bim= (
  f"{prscs_bim_prefix}.bim"
  )

prscs_bim_hm3_prefix = (
  f"{output_dir}/{cohort}_prscs_hm3")

prscs_bim_hm3= (
  f"{prscs_bim_hm3_prefix}.bim"
  )


############################################################
#Define which genotype format and which command to use
############################################################

if genotype_format == "plink1" and genotype_layout == "chromosome_split":

  make_bim_command = r"""

    echo "Creating PRS-CS BIM from chromosome-split PLINK1 files"

    cat {input.variants} > {output.bim}

  """

elif genotype_format == "plink1" and genotype_layout == "genomewide":

  make_bim_command = r"""

    echo "Using genome-wide PLINK1 BIM for PRS-CS"

    cp {input.variants} {output.bim}

  """

elif genotype_format == "plink2" and genotype_layout == "chromosome_split":

  make_bim_command = r"""

    echo "Creating PRS-CS BIM from chromosome-split PLINK2 files"

    :> {params.tmpbim}

    for pvar in {input.variants}; do

      awk 'BEGIN{{OFS="\t"}}
        !/^#/ {{
          print $1, $3, 0, $2, $4, $5
          }}' "$pvar" >> {params.tmpbim}

    done

    mv {params.tmpbim} {output.bim}
  """

elif genotype_format == "plink2" and genotype_layout == "genomewide":

  make_bim_command = r"""

    echo "Creating PRS-CS BIM from genome-wide PLINK2 PVAR"

    awk 'BEGIN{{OFS="\t"}}
      !/^#/ {{
        print $1, $3, 0, $2, $4, $5
        }}' {input.variants} > {output.bim}
  """

else: 
  raise ValueError(
    r"Unsupported genotype combination: ",
    f"{genotype_format}, {genotype_layout}"
  )

############################################################
#Define the final outputs
############################################################

final_outputs = [
  f"{output_dir}/PRS_metabolite_associations_{cohort}_{covariates_name}.txt"
  ]

if exclude_region:
  final_outputs.append(
    f"{output_dir}/PRS_metabolite_associations_{cohort}_{exclusion_suffix}_{covariates_name}.txt"
    )

rule all:
  input:
    final_outputs


############################################################
#Subet GWAS to relevant columns 
############################################################


rule extract_gwas_columns:
    input:
        sst=lambda wc: get_run(wc)["sst_file"]
    output:
        sst=f"{gwas_dir}/{{phenotype}}.tsv"
    params:
        snp_col=lambda wc: get_run(wc)["snp_col"],
        effect_allele_col=lambda wc: get_run(wc)["effect_allele_col"],
        other_allele_col=lambda wc: get_run(wc)["other_allele_col"],
        effect_col=lambda wc: get_run(wc)["effect_col"],
        effect_type=lambda wc: get_run(wc)["effect_type"],
        stat_col=lambda wc: get_run(wc)["stat_col"],
        stat_type=lambda wc: get_run(wc)["stat_type"]
    conda:
        "envs/r_base.yaml"
    script:
        "scripts/format_GWAS.R"


############################################################
#Create one combined BIM for PRS-CS
############################################################

rule make_bim: 
  input:
    variants=get_variant_files()
  output:
    bim=prscs_bim
  params:
    outdir=output_dir,
    tmpdir=f"{output_dir}/tmp_{cohort}",
    tmpbim=f"{output_dir}/tmp_{cohort}/tmp.bim",
    tmp_prefix=f"{output_dir}/tmp_{cohort}/genomewide"
  conda:
    "envs/plink.yaml"
  shell: r"""

    set -euo pipefail

    mkdir -p {params.outdir}
    mkdir -p {params.tmpdir}

  
    """ + make_bim_command + r"""
    echo "PRS-CS BIM created: {output.bim}"
    """

############################################################
#Align BIM Alleles to LD_reference to ensure cross-cohort harmonisation
############################################################

rule subset_bim_hm3:
  input:
    bim=rules.make_bim.output.bim,
    hapmap=hapmap_snps
  output:
    bim=prscs_bim_hm3,
    mapped=f"{output_dir}/{cohort}_hm3_map.txt" 
  conda:
    "envs/r_base.yaml"
  script:
    "scripts/harmonise_bim_hm3.R"

############################################################
#Run PRSCS - here one run per phenotype/prscs settings
############################################################
rule run_prscs:
  input:
    ld_ref=ld_ref,
    bim=rules.subset_bim_hm3.output.bim, 
    gwas=rules.extract_gwas_columns.output.sst
  output:
    log=f"{output_dir}/{cohort}_{{phenotype}}_phi_auto/{cohort}_{{phenotype}}_phi_auto_pst_eff_a1_b0.5_phiauto_chr{{chr}}.log",
    weights=f"{output_dir}/{cohort}_{{phenotype}}_phi_auto/{cohort}_{{phenotype}}_phi_auto_pst_eff_a1_b0.5_phiauto_chr{{chr}}.txt",
  params:
    bim=prscs_bim_hm3_prefix,
    prscs=f"PRScs/PRScs.py",
    n_gwas=lambda wc: get_run(wc)["n_gwas"],
    outdir=lambda wc: f"{output_dir}/{cohort}_{wc.phenotype}_phi_auto",
    prefix=lambda wc: f"{output_dir}/{cohort}_{wc.phenotype}_phi_auto/{cohort}_{wc.phenotype}_phi_auto"
  conda:
    "envs/prscs.yaml"
  shell: r"""
    
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
#----------------------------------------------------------------------
# BRANCH1 - for chromosome split files 
#----------------------------------------------------------------------

rule score_PRS_chr:
  input:
    genotype=get_genotype_files,
    weights=rules.run_prscs.output.weights
  output:
    sscore=f"{output_dir}/{cohort}_{{phenotype}}_phi_auto/{cohort}_{{phenotype}}_phi_auto_chr{{chr}}.sscore"
  params:
    score_dir=lambda wc: (
      f"{output_dir}/{cohort}_{wc.phenotype}_phi_auto"
      ),
    genotype_prefix=lambda wc: get_genotype_prefix(wc.chr),
    input_flag=plink_input_flag,
    outprefix=lambda wc: (
      f"{output_dir}/{cohort}_{wc.phenotype}_phi_auto/{cohort}_{wc.phenotype}_phi_auto_chr{wc.chr}"
      )
  conda:
    "envs/plink.yaml"
  shell: r"""
    set -euo pipefail
    
    mkdir -p {params.score_dir}

    echo "cohort:     {cohort}"
    echo "weights:    {input.weights}"
    echo "phenotype:  {wildcards.phenotype}"

    plink2 \
      {params.input_flag} {params.genotype_prefix} \
        --score {input.weights} 2 4 6 no-mean-imputation cols=+denom,+scoresums \
        --out {params.outprefix}

    echo  "calculating PRS for {wildcards.phenotype} for chr{wildcards.chr}"
    """

############################################################
# if regional exclusion is on, also generate scores with region excluded
############################################################

if exclude_region:
  rule score_PRS_chr_excl:
    input:
      genotype=get_genotype_files,
      weights=rules.run_prscs.output.weights,
      exclude=exclude_variants
    output:
      sscore=f"{output_dir}/{cohort}_{{phenotype}}_phi_auto/{cohort}_{{phenotype}}_phi_auto_chr{{chr}}_{exclusion_suffix}.sscore"
    params:
      score_dir=lambda wc: (
        f"{output_dir}/{cohort}_{wc.phenotype}_phi_auto"
        ),
      genotype_prefix=lambda wc: get_genotype_prefix(wc.chr),
      input_flag=plink_input_flag,
      outprefix=lambda wc: (
        f"{output_dir}/{cohort}_{wc.phenotype}_phi_auto/{cohort}_{wc.phenotype}_phi_auto_chr{wc.chr}_{exclusion_suffix}"
        )
    conda:
      "envs/plink.yaml"
    shell: r"""
      set -euo pipefail
      
      mkdir -p {params.score_dir}

      echo "cohort:     {cohort}"
      echo "weights:    {input.weights}"
      echo "phenotype:  {wildcards.phenotype}"
      echo "exclude: {exclude_locus}"

      plink2 \
        {params.input_flag} {params.genotype_prefix} \
          --exclude {input.exclude} \
          --score {input.weights} 2 4 6 no-mean-imputation cols=+denom,+scoresums \
          --out {params.outprefix}

      echo  "calculating PRS excluding {exclude_locus} for {wildcards.phenotype} for chr{wildcards.chr}"
      """

#----------------------------------------------------------------------
# BRANCH2 - for genome-wide files
#----------------------------------------------------------------------

rule combine_prscs_weights:
  input:
    weights=lambda wc: expand(
      rules.run_prscs.output.weights,
      phenotype=wc.phenotype,
      chr=chromosomes
      )
  output:
    weights=f"{output_dir}/{cohort}_{{phenotype}}_phi_auto/{cohort}_{{phenotype}}_phi_auto_weights.txt"
  shell: r"""
    set -euo pipefail

    cat {input.weights} > {output.weights}

    echo "Combined PRS-CS weights for {wildcards.phenotype}"
    """

rule score_PRS_genomewide:
  input:
    genotype=get_genomewide_genotype_files,
    weights=rules.combine_prscs_weights.output.weights
  output:
    sscore=f"{output_dir}/{cohort}_{{phenotype}}_phi_auto/{cohort}_{{phenotype}}_phi_auto.sscore"
  params:
    score_dir=lambda wc: (f"{output_dir}/{cohort}_{wc.phenotype}_phi_auto"),
    genotype_prefix=genotype_prefix,
    input_flag=plink_input_flag,
    outprefix=lambda wc: (f"{output_dir}/{cohort}_{wc.phenotype}_phi_auto/{cohort}_{wc.phenotype}_phi_auto")
  conda:
    "envs/plink.yaml"
  shell: r"""
    set -euo pipefail

    mkdir -p {params.score_dir}

    echo "cohort:     {cohort}"
    echo "weights:    {input.weights}"
    echo "phenotypes:   {wildcards.phenotype}"

    plink2 \
      {params.input_flag} {params.genotype_prefix} \
      --score {input.weights} 2 4 6 no-mean-imputation cols=+denom,+scoresums \
      --out {params.outprefix}

    echo "Calculated genome-wide PRS for {wildcards.phenotype}"
    """

############################################################
# if regional exclusion is on, also generate scores with region excluded
############################################################

if exclude_region:
  rule score_PRS_genomewide_excl:
    input:
      genotype=get_genomewide_genotype_files,
      weights=rules.combine_prscs_weights.output.weights,
      exclude=exclude_variants
    output:
      sscore=f"{output_dir}/{cohort}_{{phenotype}}_phi_auto/{cohort}_{{phenotype}}_phi_auto_{exclusion_suffix}.sscore"
    params:
      score_dir=lambda wc: (f"{output_dir}/{cohort}_{wc.phenotype}_phi_auto"),
      genotype_prefix=genotype_prefix,
      input_flag=plink_input_flag,
      outprefix=lambda wc: (f"{output_dir}/{cohort}_{wc.phenotype}_phi_auto/{cohort}_{wc.phenotype}_phi_auto_{exclusion_suffix}")
    conda:
      "envs/plink.yaml"
    shell: r"""
      set -euo pipefail

      mkdir -p {params.score_dir}

      echo "cohort:     {cohort}"
      echo "weights:    {input.weights}"
      echo "phenotypes:   {wildcards.phenotype}"
      echo "exclude: {exclude_locus}"

      plink2 \
        {params.input_flag} {params.genotype_prefix} \
        --exclude {input.exclude} \
        --score {input.weights} 2 4 6 no-mean-imputation cols=+denom,+scoresums \
        --out {params.outprefix}

      echo "Calculated genome-wide PRS excluding {exclude_locus} for {wildcards.phenotype}"
      """

############################################################
#Summarise scores - this is necessary if scores were calculated per chromosome
############################################################

rule sum_prscs:
  input:
    sscore=get_score_files
  output:
    profile=f"{output_dir}/{cohort}_{{phenotype}}_phi_auto/{cohort}_{{phenotype}}_phi_auto.profile"
  conda:
    "envs/r_base.yaml"
  script:
    "scripts/sum_scores_p12.R"



############################################################
# if regional exclusion is on, also generate scores with region excluded
############################################################

if exclude_region:
  rule sum_prscs_excl:
    input:
      sscore=get_excl_score_files
    output:
      profile=f"{output_dir}/{cohort}_{{phenotype}}_phi_auto/{cohort}_{{phenotype}}_phi_auto_{exclusion_suffix}.profile"
    conda:
      "envs/r_base.yaml"
    script:
      "scripts/sum_scores_p12.R"

############################################################
#Run association models
############################################################

rule association_models: 
  input:
    profiles=expand(
      rules.sum_prscs.output.profile,
      phenotype=df["phenotype"].tolist()
      ),
    metabolites=association["metabolites"],
    metadata=association["metadata"]
  output:
    results=f"{output_dir}/PRS_metabolite_associations_{cohort}_{covariates_name}.txt"
  params:
    covar=association["covariates"],
    covar_name=association["covariate_name"],
    phenotypes=df["phenotype"].tolist(),
  conda:
    "envs/r_base.yaml"
  script:
    "scripts/prs_metabolite_associations_p12.R"

############################################################
# if regional exclusion is on, also generate scores with region excluded
############################################################

if exclude_region:
  rule association_models_excl: 
    input:
      profiles=expand(
        rules.sum_prscs_excl.output.profile,
        phenotype=df["phenotype"].tolist()
        ),
      metabolites=association["metabolites"],
      metadata=association["metadata"]
    output:
      results=f"{output_dir}/PRS_metabolite_associations_{cohort}_{exclusion_suffix}_{covariates_name}.txt"
    params:
      covar=association["covariates"],
      covar_name=association["covariate_name"],
      phenotypes=df["phenotype"].tolist(),
    conda:
      "envs/r_base.yaml"
    script:
      "scripts/prs_metabolite_associations_p12.R"
