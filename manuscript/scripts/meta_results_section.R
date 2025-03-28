# Fill in some missing gaps needed to write up the meta-analysis results section

library(tidyr)
library(data.table)
library(glue)
library(dplyr)
library(gwascat)
library(ieugwasr)

# ------------------------------
# Get numbers for significant SNPs and numbers of genomic regions these are in.
# Paste as text to go in manuscript.

print_clumps_sentences_mrmega <- function(file_name){
  clumps <- fread(here::here("meta", "antidep-2501", file_name))
  # Only get significant SNPs
  names(clumps)[names(clumps) == "P-value_association"] <- "P"
  
  clumps <- clumps %>%
    filter(P <= 5e-8) 
  
  nSNPs <- nrow(clumps)
  nGenomicRegions <- length(unique(clumps$Locus))
  glue("{nSNPs} significant independent SNPs in {nGenomicRegions} genomic regions")
}

### MR-MEGA
lapply(c("antidep-2501-mrmega-N06A-DIV.clumps.tsv",
         "antidep-2501-mrmega-N06AA-DIV.clumps.tsv",
         "antidep-2501-mrmega-N06AB-DIV.clumps.tsv"),
       print_clumps_sentences_mrmega)

#### Fixed

get_clumps <- function(file_name){
  clumps <- fread(here::here("meta", "antidep-2501", file_name))
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

lapply(c("antidep-2501-fixed-N06A-AFR.clumps.tsv",
          "antidep-2501-fixed-N06A-EAS.clumps.tsv",
         "antidep-2501-fixed-N06A-EUR.clumps.tsv",
         "antidep-2501-fixed-N06A-SAS.clumps.tsv",
         "antidep-2501-fixed-N06AA-AFR.clumps.tsv",
         "antidep-2501-fixed-N06AA-EUR.clumps.tsv",
         "antidep-2501-fixed-N06AB-AFR.clumps.tsv",
         "antidep-2501-fixed-N06AB-EUR.clumps.tsv"),
       print_clumps_sentences_fixed)

# ------------------------------
# Run the GWAS catalogue look up code for the fixed effects results
gwcat <- get_cached_gwascat()

clumps <- lapply(c("antidep-2501-fixed-N06A-AFR.clumps.tsv",
                   "antidep-2501-fixed-N06A-EAS.clumps.tsv",
                   "antidep-2501-fixed-N06A-EUR.clumps.tsv",
                   "antidep-2501-fixed-N06A-SAS.clumps.tsv",
                   "antidep-2501-fixed-N06AA-AFR.clumps.tsv",
                   "antidep-2501-fixed-N06AA-EUR.clumps.tsv",
                   "antidep-2501-fixed-N06AB-AFR.clumps.tsv",
                   "antidep-2501-fixed-N06AB-EUR.clumps.tsv"),
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

names(gwascat_tables) <- c("N06A-AFR", "N06A-EAS", "N06A-EUR", "N06A-SAS", "N06AA-AFR", "N06AA-EUR", "N06AB-AFR", "N06AB-EUR")

save_gwascat_table <- function(gwascat_table, file_name){
  write.csv(gwascat_table, 
            here::here("manuscript", "tables", paste0("gwascat_fixed_table_", file_name, ".csv")), 
            quote = F, row.names = F)
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



