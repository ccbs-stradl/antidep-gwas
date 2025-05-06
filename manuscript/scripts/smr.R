# Create .cols sidecar meta data for SMR results

# load libraries
library(dplyr)
library(tidyr)
library(stringr)

# Create col name descriptions, same for all SMR results
colname_descriptions <-
c("Gene" = "placeholder",
"qtl_name" = "placeholder",
"probeID" = "placeholder",
"ProbeChr" = "placeholder",
"Probe_bp" = "placeholder",
"topSNP" = "placeholder",
"topSNP_chr" = "placeholder",
"topSNP_bp" = "placeholder",
"A1" = "placeholder",
"A2" = "placeholder",
"Freq" = "placeholder",
"b_GWAS" = "placeholder",
"se_GWAS" = "placeholder",
"p_GWAS" = "placeholder",
"b_eQTL" = "placeholder",
"se_eQTL" = "placeholder",
"p_eQTL" = "placeholder",
"b_SMR" = "placeholder",
"se_SMR" = "placeholder",
"p_SMR" = "placeholder",
"p_HEIDI" = "placeholder",
"nsnp_HEIDI" = "placeholder",
"gene_id" = "placeholder",
"chr" = "placeholder",
"start" = "placeholder",
"end" = "placeholder",
"strand" = "placeholder",
"GWAS_LOCUS" = "placeholder",
"Lead_SNP" = "placeholder",
"Lead_SNP_BP" = "placeholder"
)

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
