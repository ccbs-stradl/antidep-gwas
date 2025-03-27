# Fill in some missing gaps needed to write up the meta-analysis results section

library(data.table)
library(glue)
library(dplyr)

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
print_clumps_sentences_fixed <- function(file_name){
  clumps <- fread(here::here("meta", "antidep-2501", file_name))
  # Only get significant SNPs
  clumps <- clumps %>%
    filter(P <= 5e-8) 
  
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

gwascat_N06A_table <- look_up_snps(clumps_N06A, gwcat)
gwascat_N06AA_table <- look_up_snps(clumps_N06AA, gwcat)
gwascat_N06AB_table <- look_up_snps(clumps_N06AB, gwcat)

write.csv(gwascat_N06A_table, "manuscript/tables/gwas_cat_N06A_table.csv", quote = F, row.names = F)

# ------------------------------
# Get the number of unique SNPs and check for duplicate trait names in GWAS catalogue look up

