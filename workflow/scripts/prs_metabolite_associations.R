###############################################################################
# PRS-lipid association testing
# Integrated into the PRScs Snakemake pipeline.
# Author: jaenice94  
# Date: 2024-06-15
# Runs once as an aggregation step: reads every scored .profile, builds the
# PRS matrix, and fits lipid ~ PRS + covariates & lipid ~ covariates models
###############################################################################

library(data.table);library(dplyr)

profiles_vec <- unlist(snakemake@input[["profiles"]]) #read .profiles
phenotypes <- unlist(snakemake@params[["phenotypes"]]) #list of phenotypes
metabolites_file <- snakemake@input[["metabolites"]]
metadata_file <- snakemake@input[["metadata"]]
covariates <- snakemake@params[["covar"]]

#order is as in samplesheet 
names(profiles_vec) <- phenotypes

metabolites <- fread(metabolites_file)
metadata <- fread(metadata_file)

metabolite_meta <- merge(metabolites, metadata, by="IID")

#merge PRS-profiles into one file
PRS <- NULL
for (p in phenotypes) {
	dt <- fread(profiles_vec[[p]]) %>%
		transmute(IID, !!p := SCORESUM)
	PRS <- if (is.null(PRS)) dt else inner_join(PRS, dt, by = "IID")
}

metabolite_meta_PRS <- merge(metabolite_meta, PRS, by = "IID")

prs_map <- setNames(phenotypes, phenotypes)

lipid_names <- setdiff(colnames(metabolites), c("IID", "FID")) #if no FID, will ignore

# baseline:  metabolite ~ covariates
fit_baseline <- function(lipid_names, data, covar_terms) {
	rhs <- paste(covar_terms, collapse = " + ") #extract covariates and bring to formula format

	#for each lipid build a lm with this formula; combine all the data.frames into one
	do.call(rbind, lapply(lipid_names, function(x) {
		s <- summary(lm(as.formula(paste(x, "~", rhs)), data = data)) 
		data.frame(model=x, Adj_r_sq_baseline = s$adj.r.squared, stringsAsFactors = FALSE) 
	}))
}

# prs model:  metabolite ~ prs + covariates
fit_prs <- function(prs_col, lipid_names, data, covar_terms){
	rhs <- paste(c(prs_col,covar_terms), collapse = " + ") #extract covariates +prs and bring to formula format
	do.call(rbind, lapply(lipid_names, function(x) {
		s <- summary(lm(as.formula(paste(x, "~", rhs)), data = data)) 
		co <- s$coefficients
		data.frame(model=x,
				Effect_PRS = co[prs_col, "Estimate"],
				P_value_PRS = co[prs_col, "Pr(>|t|)"],
				Adj_r_sq_model=s$adj.r.squared,
				stringsAsFactors = FALSE)
	}))
}

#function to run both and calculate r2
run_prs_lm <- function(data, lipid_names, prs_map, covars, covar_labels){

	baseline <- fit_baseline(lipid_names, data, covars)

	run_one_prs <- function(label, prs_col){
		fit_prs(prs_col, lipid_names, data, covars) %>%
		merge(baseline, by = "model") %>%
		mutate(r_squared_diff = Adj_r_sq_model - Adj_r_sq_baseline, 
				PRS = label)
	}

	#combine outputs 

	bind_rows(lapply(names(prs_map), function(lbl) run_one_prs(lbl, prs_map[[lbl]]))) %>%
		mutate(logP = -log10(P_value_PRS))
}

covars_base <- covariates

results <- run_prs_lm(
	data=metabolite_meta_PRS,
	lipid_names=lipid_names,
	prs_map=prs_map,
	covars=covars_base
	)

fwrite(results, snakemake@output[["results"]], sep = "\t")
