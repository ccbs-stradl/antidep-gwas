# Fill in some missing gaps needed to write up the meta-analysis results section

library(tidyr)
library(data.table)
library(glue)
library(dplyr)
library(gwascat)
library(ieugwasr)
library(stringr)
library(readr)

# ------------------------------
# Get numbers for significant SNPs and numbers of genomic regions these are in.
# Paste as text to go in manuscript.

print_clumps_sentences_mrmega <- function(file_name){
  clumps <- fread(file_name)
  # Only get significant SNPs
  names(clumps)[names(clumps) == "P-value_association"] <- "P"
  
  clumps <- clumps %>%
    filter(P <= 5e-8) 
  
  nSNPs <- nrow(clumps)
  nGenomicRegions <- length(unique(clumps$Locus))
  glue("{nSNPs} significant independent SNPs in {nGenomicRegions} genomic regions")
}

### MR-MEGA
path <- here::here("results", "meta", "antidep-2501")
files <- list.files(path, ".+mrmega.+clumps\\.tsv", full.names = T)
lapply(files,
       print_clumps_sentences_mrmega)

#### Fixed

get_clumps <- function(file_name){
  clumps <- fread(file_name)
  # Only get significant SNPs
  clumps <- clumps %>%
    filter(P <= 5e-8) 
}

print_clumps_sentences_fixed <- function(file_name){
  clumps <- get_clumps(file_name)
  
  nSNPs <- nrow(clumps)
  nGenomicRegions <- length(unique(clumps$LOCUS))
  glue("{nSNPs} significant independent SNPs in {nGenomicRegions} genomic regions")
}

path <- here::here("results", "meta", "antidep-2501")
files <- list.files(path, ".+fixed.+clumps\\.tsv", full.names = T)
lapply(files,
       print_clumps_sentences_fixed)

# ------------------------------
# Run the GWAS catalogue look up code for the fixed effects results
gwcat <- get_cached_gwascat()

# if the above doesn't work, try:
# library(BiocFileCache)
# 
# # Create a BiocFileCache instance
# bfc <- BiocFileCache::BiocFileCache()
# 
# # Add the GWAS catalog URL manually to the cache
# bfcadd(
#   bfc,
#   rname = "gwascat",
#   fpath = "http://www.ebi.ac.uk/gwas/api/search/downloads/alternative",
#   rtype = "web"
# )

clumps <- lapply(files,
                 get_clumps)

look_up_snps <- function(clump, gwcat){
  open_gwas <- phewas(variants = clump |> as_tibble() |> pull(ID), pval=5e-8)
  
  gwcat_snps <-
    gwcat |>
    select(PUBMEDID, `DISEASE/TRAIT`, SNPS, MAPPED_TRAIT) |>
    separate_wider_delim(SNPS, delim = "; ", names_sep = "_", too_few = "align_start") |>
    pivot_longer(starts_with("SNP"), values_to = 'SNP') |>
    filter(!is.na(SNP)) |>
    select(-name)
  
  nomhc_snps <- clump |> as_tibble() |> pull(ID)
  
  gwcat_snps |>
    filter(SNP %in% nomhc_snps) |>
    transmute(`DISEASE/TRAIT` = str_sub(`DISEASE/TRAIT`, 1, 50), SNP) |>
    arrange(SNP)
}

gwascat_tables <- lapply(clumps, look_up_snps, gwcat = gwcat)

# Use regular expression to rename the tables
# "[A-Za-z0-9]+" is a regex pattern that matches one or more alphanumeric characters
# "-" matches the literal hyphen between the two alphanumeric strings
# "(?=)" is a lookahead that ensures the preceding match is 
# followed by the given pattern but does not include it in the match.
# "." are escapped wtih "\\" to match the literal period in the file name
names(gwascat_tables) <- str_extract(files, "([A-Za-z0-9]+-[A-Za-z0-9]+)(?=\\.clumps\\.tsv)")

save_gwascat_table <- function(gwascat_table, file_name){
  gwascat_colname_descriptions <- 
    c("DISEASE/TRAIT" = "name of disease or trait phenotype",
      "SNP" = "variant identifier"
    )
  
  colname_descriptions_table <- tibble(column = names(gwascat_colname_descriptions), description = gwascat_colname_descriptions)
  
  write.csv(gwascat_table, 
            here::here("manuscript", "tables", paste0("gwascat_fixed_table_", file_name, ".csv")), 
            quote = F, row.names = F)
  
  if(any(colname_descriptions_table$column != colnames(gwascat_table))){
    stop(glue("Column names in {file_name_cs} are not all described in colname_descriptions"))
  }
  
  write_tsv(colname_descriptions_table,here::here("manuscript", "tables", paste0("gwascat_fixed_table_", file_name, ".csv.cols")))
}

Map(save_gwascat_table, gwascat_tables, names(gwascat_tables))

# ------------------------------
# Get the number of unique SNPs and check for duplicate trait names in GWAS catalogue look up

# fixed
lapply(gwascat_tables, function(x) {
  nSNPs <- x %>% pull(SNP) %>% unique() %>% length()
  nDiseases <- x %>% pull(`DISEASE/TRAIT`) %>% unique() %>% length()
  glue("{nSNPs} unique SNPs were identified in the NHGRI-EBI GWAS catalogue as associated with {nDiseases} traits ")
})

# MR-MEGA
gwas_cat_results <- fread(here::here("scripts", "multi_files", "gwas_cat_N06A_table.csv"))
nSNPs <- gwas_cat_results %>% pull(SNP) %>% unique() %>% length()
nDiseases <- gwas_cat_results %>% pull(`DISEASE/TRAIT`) %>% unique() %>% length()
glue("{nSNPs} unique SNPs were identified in the NHGRI-EBI GWAS catalogue as associated with {nDiseases} traits ")

gwas_cat_results <- fread(here::here("scripts", "multi_files", "gwas_cat_N06AA_table.csv"))
nSNPs <- gwas_cat_results %>% pull(SNP) %>% unique() %>% length()
nDiseases <- gwas_cat_results %>% pull(`DISEASE/TRAIT`) %>% unique() %>% length()
glue("{nSNPs} unique SNPs were identified in the NHGRI-EBI GWAS catalogue as associated with {nDiseases} traits ")

gwas_cat_results <- fread(here::here("scripts", "multi_files", "gwas_cat_N06AB_table.csv"))
nSNPs <- gwas_cat_results %>% pull(SNP) %>% unique() %>% length()
nDiseases <- gwas_cat_results %>% pull(`DISEASE/TRAIT`) %>% unique() %>% length()
glue("{nSNPs} unique SNPs were identified in the NHGRI-EBI GWAS catalogue as associated with {nDiseases} traits ")



