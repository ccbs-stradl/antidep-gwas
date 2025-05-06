# Create .cols sidecar meta data for SMR results

# load libraries
library(dplyr)
library(tidyr)
library(stringr)

# Create col name descriptions, same for all SMR results
# colname_descriptions <- 
# c(index
# qtl_name
# probeID
# ProbeChr
# Probe_bp
# topSNP
# topSNP_chr
# topSNP_bp
# A1
# A2
# Freq
# b_GWAS
# se_GWAS
# p_GWAS
# b_eQTL
# se_eQTL
# p_eQTL
# b_SMR
# se_SMR
# p_SMR
# p_HEIDI
# nsnp_HEIDI
# gene_id
# chr
# start
# end
# strand
# GWAS_LOCUS
# Lead_SNP
# Lead_SNP_BP)

# Load function to create .cols sidecar meta data file
source(here::here("manuscript/scripts/supplementary_tables_excell_create_cols_meta_FUN.R"))

# Read in all SMR results for blood and brain
paths <- list(
  "results/maps/smr/blood/trait_eSMR.merged.tsv",
  "results/maps/smr/blood/trait_mSMR.merged.tsv",
  "results/maps/smr/blood/trait_pSMR.merged.tsv",
  "results/maps/smr/brainmeta/trait_eSMR.merged.tsv",
  "results/maps/smr/brainmeta/trait_mSMR.merged.tsv",
  "results/maps/smr/brainmeta/trait_sSMR.merged.tsv"
)

lapply(paths, function(path){
  smr_results <- fread(here::here(path))
  
  create_cols_meta(
    file_name = path,
    table_variable_name = smr_results,
    colname_descriptions = colname_descriptions
  )
})
