# GWAS catalog look up for MDD (Adams et al. 2025)

library(tidyr)
library(data.table)
library(glue)
library(dplyr)
library(gwascat)
library(purrr)
library(ieugwasr)
library(stringr)
library(readr)
library(readxl)

# Download COJO results from: https://figshare.com/articles/online_resource/Supplementary_Results_for_GWAS_meta-analysis_of_major_depression_PGC_MDD2025_/27089614?file=54826145

################################
######## FUNCTIONS #############
read_sheet <- function(sheet_index, path) {
  read_xlsx(here::here('manuscript/MDD_GWAS_Online_Results_(COJO).xlsx'), sheet = sheet_index)
}

pull_SNPs <- function(sheet, p_threshold = 5e-8){
  sheet |>
    filter(P <= p_threshold) |>
    pull(SNP)
}

get_MDD_SNPs <- function(sheet_index, path){
  sheet <- read_sheet(sheet_index, path)
  SNPs <- pull_SNPs(sheet)
}

look_up_snps <- function(SNPs, gwcat){
  gwcat_snps <-
    gwcat |>
    select(PUBMEDID, `DISEASE/TRAIT`, SNPS, MAPPED_TRAIT) |>
    separate_wider_delim(SNPS, delim = "; ", names_sep = "_", too_few = "align_start") |>
    pivot_longer(starts_with("SNP"), values_to = 'SNP') |>
    filter(!is.na(SNP)) |>
    select(-name)
  
  gwcat_snps |>
    filter(SNP %in% SNPs) |>
    transmute(`DISEASE/TRAIT` = str_sub(`DISEASE/TRAIT`, 1, 50), SNP) |>
    arrange(SNP)
}

main <- function(gwcat, ancestry, path){
  SNPs <- lapply(ancestry, get_MDD_SNPs, path = path)
  gwascat_tables <- lapply(SNPs, look_up_snps, gwcat = gwcat)
  gwascat_tables
}

################################
######## Global variables ######
if(!exists("gwcat")) {
  gwcat <- get_cached_gwascat()
}

# Name the sheet index with the ancestry
# (using a named list here ensures the SNPs variable is named later)
ancestry <- list("multi" = 2,
                 "EUR" = 3)

gwascat_tables  <- main(gwcat = gwcat, ancestry = ancestry, path = here::here('manuscript/MDD_GWAS_Online_Results_(COJO).xlsx'))
gwascat_tables

# ------------------------------
