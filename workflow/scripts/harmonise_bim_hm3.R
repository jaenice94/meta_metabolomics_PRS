###############################################################################
# Harmonise bim to hapmap snps 
# Integrated into the PRScs Snakemake pipeline.
# Author: jaenice94  
# Date: 2026-09-02
#
# PRS-CS correctly handles allele swaps and strand orientation when matching
# variants. However, changing allele orientation changes the sign of the
# corresponding effect representation. Because posterior effects are sampled
# stochastically during MCMC, using the same random seed after an allele
# sign reversal does not produce an exactly sign-reversed MCMC trajectory.
#
# Consequently, equivalent genotype datasets represented with different
# A1/A2 orientations can yield different MCMC trajectories and estimates of
# the global shrinkage parameter.
#
# To maximise reproducibility and comparability across cohorts, cohort BIM
# files are therefore restricted and oriented consistently to the allele
# representation in the PRS-CS HapMap3 LD reference.
#
# Variants are classified into:
# same, swap, strand, strand_swap, ambiguous, mismatch
###############################################################################

library(data.table);library(dplyr)

bim_file <- snakemake@input[["bim"]]
hapmap_file <- snakemake@input[["hapmap"]]
output_file <- snakemake@output[["bim"]]
output_map <- snakemake@output[["mapped"]]

bim <- fread(bim_file, 
             header = FALSE,
             col.names = c(
               "CHR",
               "SNP",
               "CM",
               "BP",
               "A1",
               "A2"
             )
            ) %>%
  mutate(
    BIM_ORDER = row_number()
  )

hapmap <- fread(hapmap_file, header = TRUE) %>%
  transmute(
    SNP,
    REF_CHR = CHR,
    REF_BP = BP,
    REF_A1 = A1,
    REF_A2 = A2
  )

#subset bim by rsid 
bim_hapmap <- bim %>% filter(SNP %in% hapmap$SNP)

#report:
message(
  "Original BIM variants: ",
  nrow(bim)
)

if (anyDuplicated(bim$SNP)) {
  stop(
    "Duplicated SNP IDs found in cohort BIM. Please filter out duplicate IDs!"
  )
}

message(
  "Variants in HapMap3: ",
  nrow(hapmap)
)

message(
  "Variants in common between BIM and HapMap3: ",
  nrow(bim_hapmap)
)


################## Harmonise

harm <- bim_hapmap %>% inner_join(hapmap, by = "SNP")

harm <- harm %>%
  mutate(
    REF_A1_COMP = recode(
      REF_A1,
      A = "T",
      T = "A",
      C = "G",
      G = "C"
    ),
    REF_A2_COMP = recode(
      REF_A2,
      A = "T",
      T = "A",
      C = "G",
      G = "C"
    )
  )

valid_alleles <- c("A", "C", "G", "T")

harm <- harm %>%
  mutate(
    VALID_SNP = 
      A1 %in% valid_alleles &
      A2 %in% valid_alleles &
      REF_A1 %in% valid_alleles &
      REF_A2 %in% valid_alleles,
    
    AMBIGUOUS =
      paste0(REF_A1, REF_A2) %in% c("AT", "TA", "CG", "GC"),
    
    CHR_MATCH =
      as.character(CHR) == as.character(REF_CHR)
  )

#Classify

harm <- harm %>% 
  mutate(
    STATUS = case_when(
      !CHR_MATCH ~ "chr_mismatch",
      !VALID_SNP ~ "non_acgt",
      AMBIGUOUS ~ "ambiguous",
      A1 == REF_A1 & A2 == REF_A2 ~ "same",
      A1 == REF_A2 & A2 == REF_A1 ~ "swap",
      A1 == REF_A1_COMP & A2 == REF_A2_COMP ~ "strand",
      A1 == REF_A2_COMP & A2 == REF_A1_COMP ~ "strand_swap",
      TRUE ~ "mismatch"
    )
)

message("Harmonisation status: ")

status_counts <- harm %>% count(STATUS, name = "N") %>% arrange(STATUS)

print(status_counts)

keep_status <- c("same", "swap", "strand", "strand_swap")

harm_keep <- harm %>% filter(STATUS %in% keep_status)

message(
  "Variants retained after harmonisation: ",
  nrow(harm_keep))

message(
  "Variants dropped after harmonisation: ",
  nrow(harm)-nrow(harm_keep))

# for the final bim - either use REF_A1 / REF_A2; or if strand flipped, keep original allele combination, but orient in respect to hapmap -> avoid scoring issue

harm_keep <- harm_keep %>%
  mutate(
    ALIGNED_A1 = case_when(
      STATUS %in% c("same", "strand") ~ A1,
      STATUS %in% c("swap", "strand_swap") ~ A2),
    ALIGNED_A2 = case_when(
      STATUS %in% c("same", "strand") ~ A2,
      STATUS %in% c("swap", "strand_swap") ~ A1)
  )

bim_harmonised <- harm_keep %>%
  arrange(BIM_ORDER) %>%
  transmute(
    CHR,
    SNP,
    CM,
    BP,
    A1 = ALIGNED_A1,
    A2 = ALIGNED_A2
  )

fwrite(
  bim_harmonised,
  output_file,
  sep = "\t",
  col.names = FALSE
)

harmonisation_map <- harm %>%
  arrange(BIM_ORDER) %>%
  transmute(
    SNP,
    CHR,
    BP,
    COHORT_A1 = A1,
    COHORT_A2 = A2,
    REF_CHR,
    REF_BP,
    REF_A1,
    REF_A2,
    STATUS
  )

fwrite(
  harmonisation_map,
  output_map,
  sep = "\t",
  col.names = TRUE
)
