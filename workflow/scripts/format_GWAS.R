###############################################################################
# Extract columns from GWAS for PRScs
# Integrated into the PRScs Snakemake pipeline.
# Author: jaenice94  
# Date: 2026-08-25
# based on inputs in samplesheet, extracts the necessary columns for prscs
###############################################################################

library(data.table)

# input GWAS file 
infile <- snakemake@input[["sst"]]

# output file in PRS-CS format
outfile <- snakemake@output[["sst"]]

######################################
#parameters
######################################
snp_col <- as.character(
  snakemake@params[["snp_col"]]
)

effect_allele_col <- as.character(
  snakemake@params[["effect_allele_col"]]
)

other_allele_col <- as.character(
  snakemake@params[["other_allele_col"]]
)

effect_col <- as.character(
  snakemake@params[["effect_col"]]
)

effect_type <- as.character(
  snakemake@params[["effect_type"]]
)

stat_col <- as.character(
  snakemake@params[["stat_col"]]
)

stat_type <- as.character(
  snakemake@params[["stat_type"]]
)

######################################
#Format
######################################
dt <- fread(infile)

req <- unique(c(
  snp_col,
  effect_allele_col,
  other_allele_col,
  effect_col,
  stat_col
))

missing <- setdiff(req, colnames(dt))

if (length(missing) > 0) {
  stop(
    "Missing required columns in summary statistics: ",
    paste(missing, collapse = ", ")
  )
}

######################################
#Standardise column names
######################################

prs_dt <- dt[, .(
  SNP = as.character(get(snp_col)),
  A1 = toupper(as.character(get(effect_allele_col))),
  A2 = toupper(as.character(get(other_allele_col))),
  EFFECT = get(effect_col),
  STAT = get(stat_col)
)]

setnames(
  prs_dt,
  old = c("EFFECT", "STAT"),
  new = c(effect_type, stat_type)
)

#bim creation and GWAS filtering can run in parallel, so create if doesn't exist yet
dir.create(
  dirname(outfile),
  recursive = TRUE,
  showWarnings = FALSE
)

# write out
fwrite(
  prs_dt,
  file = outfile,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

cat("Wrote PRS-CS formatted file to:\n", outfile, "\n")
