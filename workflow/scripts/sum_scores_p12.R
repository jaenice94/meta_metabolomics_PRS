###############################################################################
# Sum PRS scores across chromosomes
# Integrated into the PRScs Snakemake pipeline.
# Author: jaenice94  
# Date: 2026-08-25
#
# For chromosome-split genotype data:
#   combines chromosome-specific SCORE1_SUM values per individual.
#
# For genome-wide genotype data:
#   uses the SCORE1_SUM from the single genome-wide score file.
#
# Produces the final raw and standardized PRS used for association testing.
###############################################################################


library(data.table);library(dplyr)

files <- unlist(
  snakemake@input[["sscore"]],
  use.names = FALSE
)

output=snakemake@output[["profile"]]

if (length(files) == 0) {
  stop("No .sscore files provided.")
}

missing_files <- files[!file.exists(files)]

if (length(missing_files) > 0) {
  stop(
    "Missing .sscore files: ",
    paste(missing_files, collapse = ", ")
  )
}


##########################################
# validate and read score files
##########################################

score_tables <- lapply(files, function(f) {
  dt <- fread(f)
  
  #standardise PLINK ID column
  if ("#FID" %in% names(dt)) {
    setnames(dt, "#FID", "FID")
  }
  
  if ("#IID" %in% names(dt)) {
    setnames(dt, "#IID", "IID")
  }
  
  if (!"FID" %in% names(dt)) {
    dt[, FID := "0"]
  }
  
  #IID is required
  if (!"IID" %in% names(dt)) {
    stop("IID column missing from score files: ", f)
  }
  
  #SCORE1_SUM required for PRS calculation
  if (!"SCORE1_SUM" %in% names(dt)) {
    stop("SCORE1_SUM column missing from score files: ", f)
  }

  if (anyDuplicated(dt[, IID])) {
    stop("Duplicate individuals detected in score file: ", f)
  }
  
  dt
})


##########################################
# Combine score files
##########################################  

#if no FID column is present -> FID entries will become NA 
profiles <- score_tables %>%
  bind_rows(.id = "score_file") %>%
  group_by(FID, IID) %>%
  summarise(
    n_score_files = n_distinct(score_file),
    PRS_SUM = sum(SCORE1_SUM),
    .groups = "drop"
  )

##########################################
# Validate combined scores
##########################################

if (any(profiles$n_score_files != length(files))) {
  stop("Not all individuals are present in all score files.")
} 

if (anyNA(profiles$PRS_SUM)) {
  stop("Missing PRS values detected in score files.")
}

##########################################
# Standardise PRS 
##########################################

profiles <- profiles %>%
  mutate(
    PRS_STD = as.numeric(scale(PRS_SUM))
  ) %>%
  select(-n_score_files)
  
fwrite(profiles, output, sep = "\t")
