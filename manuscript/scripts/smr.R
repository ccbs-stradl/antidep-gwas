# Copy main results from datastore to here so they can be git tracked.
# Create .cols sidecar meta data for SMR results

# load libraries
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
library(readr)
library(glue)
library(here)

# Load function to create .cols sidecar meta data file
source(here::here("manuscript/scripts/supplementary_tables_excell_create_cols_meta_FUN.R"))

# -----------------------------------------------
# Copy results to this projects directory, from datastore
# First I made a sym link to the datastore dir (see docs/ for how to do this)
# Read in all SMR results for blood and brain
paths <- list(
  "results/maps/smr/blood/trait_eSMR.merged.tsv",
  "results/maps/smr/blood/trait_mSMR.merged.tsv",
  "results/maps/smr/blood/trait_pSMR.merged.tsv",
  "results/maps/smr/brainmeta/trait_eSMR.merged.tsv",
  "results/maps/smr/brainmeta/trait_mSMR.merged.tsv",
  "results/maps/smr/brainmeta/trait_sSMR.merged.tsv"
)

# rename file from "blood/trait_eSMR.merged.tsv" to 
# "blood_trait_eSMR.merged.tsv" etc.
# then move this file to manuscripts/tables/
rename_file_and_move <- function(old_path, new_path_prefix){
  # rename file
  tissue_type <- basename(dirname(old_path)) # extracts the last directory name
  new_file_name <- str_c(tissue_type, '_', basename(old_path))
  
  # full new path
  new_path <- paste0(new_path_prefix, "/", new_file_name)
  
  # copy results and rename with new name
  file.copy(
    from = here::here(old_path),
    to = here::here(new_path),
    overwrite = TRUE
  )
  
  return(new_path)
}

paths <- lapply(paths, rename_file_and_move, new_path_prefix = "manuscript/tables")
# -----------------------------------------------
# Anomaly - "manuscript/tables/blood_trait_pSMR.merged.tsv" should have column called Gene not index
blood_pSMR <- fread("manuscript/tables/blood_trait_pSMR.merged.tsv")
# rename
blood_pSMR <- rename(blood_pSMR, "Gene" = "index")
# rewrite
write_tsv(blood_pSMR, "manuscript/tables/blood_trait_pSMR.merged.tsv")

# -----------------------------------------------
# Create col name descriptions, same for all SMR results
colname_descriptions <-
  c("Gene" = "gene name",
    "qtl_name" = "QTL name",
    "probeID" = "probe ID",
    "ProbeChr" = "probe chromosome",
    "Probe_bp" = "probe position",
    "topSNP" = "SNP name",
    "topSNP_chr" = "SNP chromosome",
    "topSNP_bp" = "SNP position",
    "A1" = "the effect (coded) allele",
    "A2" = "the other allele",
    "Freq" = "frequency of the effect allele (estimated from the reference samples)",
    "b_GWAS" = "effect size from GWAS",
    "se_GWAS" = "standard error from GWAS",
    "p_GWAS" = "p-value from GWAS",
    "b_eQTL" = "effect size from eQTL study",
    "se_eQTL" = "standard error from eQTL study",
    "p_eQTL" = "p-value from eQTL study",
    "b_SMR" = "effect size from SMR",
    "se_SMR" = "standard error from SMR",
    "p_SMR" = "p-value from SMR",
    "p_HEIDI" = "p-value from HEIDI (HEterogeneity In Depedent Instruments) test",
    "nsnp_HEIDI" = "number of SNPs used in the HEIDI test",
    "gene_id" = "Gene ID",
    "chr" = "Chromosome",
    "start" = "Start position",
    "end" = "end position",
    "strand" = "strand",
    "GWAS_LOCUS" = "gwas locus",
    "Lead_SNP" = "lead SNP",
    "Lead_SNP_BP" = "lead SNP position"
  )

lapply(paths, function(path){
  smr_results <- fread(here::here(path))
  
  create_cols_meta(
    file_name = path,
    table_variable_name = smr_results,
    colname_descriptions = colname_descriptions
  )
})
