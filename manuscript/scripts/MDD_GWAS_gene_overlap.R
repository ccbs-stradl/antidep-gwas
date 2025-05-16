# Compare results with MDD GWAS (Adams et al. 2025)
# create a table with each row is a gene, and columns for different methods,
# with further columns for if the gene was identified in MDD GWAS (noting the method)
# ------------------------------
# Read in Supplementary Table 8 from Adams et al. 2025 (Gene mapping results summarized across methods, related to STAR Methods)
# https://www.cell.com/cell/fulltext/S0092-8674(24)01415-6

library(tidyverse)
library(tidyr)
library(glue)
library(dplyr)
library(data.table)
library(readr)
source(here::here("manuscript/scripts/supplementary_tables_excell_create_cols_meta_FUN.R"))

mdd_genes <- fread(here::here("manuscript/MDD_GWAS_supp_table_8B.csv")) %>%
  # add prefix "MDD_GWAS_" to columns: c(Nearest_gene Fine_mapping Expression Protein fastBAT HMAGMA PsyOPS)
  rename_with(~ paste0("MDD_GWAS_", .), c(Nearest_gene, Fine_mapping, Expression, Protein, fastBAT, HMAGMA, PsyOPS)) %>%
  # rename Chrom_b37 pos_start_b37 pos_end_b37 to Chr     Start       End
  rename(Chr = Chrom_b37, Start = pos_start_b37, End = pos_end_b37) %>%
  # remove cols ending in b38
  select(-ends_with("b38")) %>%
  mutate(across(c(Chr, Start, End), ~ as.character(.)))


# ------------------------------
# Read in mBAT-combo results
mBAT_combo_results <- fread(here::here("manuscript/tables/mBAT-combo.csv"))

mBAT_combo_results %>%
  group_by(phenotype, ancestry) %>%
  summarise(
    n_genes = n_distinct(Gene),
    n_genes_in_MDD = sum(Gene %in% mdd_genes$ENSID)
  )


# pivot from long table to wide table for phenotype
mBAT_combo_wide <- mBAT_combo_results %>%
  select(phenotype, Gene, gene_name, Chr, Start, End) %>%
  mutate(value = TRUE) %>% # Create a column to hold TRUE values
  pivot_wider(
    names_from = phenotype,
    values_from = value,
    names_prefix = "mBAT_combo_EUR_",
    values_fill = list(value = FALSE) # Fill missing values with FALSE
  )

# ------------------------------
# Read in susiex results
susiex_results <- data.frame(
  gene_name = c("PLCL2", "WNT3"),
  SuSiEx = TRUE
)

# Join with mbat combo results
antidep_results <- merge(mBAT_combo_wide, susiex_results, by = c("gene_name"), all = TRUE) %>%
  mutate(SuSiEx = ifelse(is.na(SuSiEx), FALSE, SuSiEx)) %>%
  mutate(across(starts_with("mBAT_combo"), ~ ifelse(is.na(.), FALSE, .)))

# Check any NA
antidep_results %>%
  filter(is.na(Gene))

# Manually fill in gene details for WNT3 (https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000108379;r=17:44839872-44910520)
antidep_results[antidep_results$gene_name == "WNT3", c("Gene", "Chr", "Start", "End")] <- c("ENSG00000108379", "17", 44839872, 44910520)

# ------------------------------
# Read in SMR results
# Relies on output from manuscript/scripts/smr.R
# this script moves and cleans SMR results to manuscripts/tables
source(here::here("manuscript/scripts/smr.R"))

paths <- list(
  "manuscript/tables/blood_trait_eSMR.merged.tsv",
  "manuscript/tables/blood_trait_mSMR.merged.tsv",
  "manuscript/tables/blood_trait_pSMR.merged.tsv",
  "manuscript/tables/brainmeta_trait_eSMR.merged.tsv",
  "manuscript/tables/brainmeta_trait_mSMR.merged.tsv",
  "manuscript/tables/brainmeta_trait_sSMR.merged.tsv"
)

smr_results_list <- lapply(paths, function(path) {
  fread(here::here(path))
})

names <- basename(unlist(paths)) %>%
  str_remove(".merged.tsv")

names(smr_results_list) <- names

# filter to get genes with p_SMR_Bonferroni < 0.05 and p_HEIDI > 0.05
smr_results <- lapply(seq_along(smr_results_list), function(i) {
  smr_results_list[[i]] %>%
    filter(p_SMR_Bonferroni < 0.05 & p_HEIDI > 0.05) %>%
    select(gene_name = Gene, Gene = gene_id, Chr = chr, Start = start, End = end) %>%
    mutate(!!sym(names[i]) := TRUE)
}) %>%
  reduce(full_join, by = c("Gene", "gene_name", "Chr", "Start", "End")) %>%
  mutate(across(ends_with("SMR"), ~ replace_na(., FALSE)))


antidep_results <- merge(antidep_results, smr_results, by = c("Gene", "gene_name", "Chr"), all = TRUE)

# Rename Start.x and Start.y to mBATcombo_SuSiEx_Start 
# and Start.y to SMR_Start
antidep_results <- antidep_results %>%
  rename(
    mBATcombo_SuSiEx_Start = Start.x,
    SMR_Start = Start.y,
    mBATcombo_SuSiEx_End = End.x,
    SMR_End = End.y
  ) %>%
  relocate(c(SMR_Start, SMR_End), .after = mBATcombo_SuSiEx_End) 

# ------------------------------
# Join with MDD GWAS results
antidep_results <- antidep_results %>%
  rename(ENSID = Gene, Gene = gene_name) %>%
  merge(mdd_genes, by = c("ENSID", "Gene", "Chr"), all = TRUE) %>%
  relocate(c(Start, End), .after = SMR_End) %>%
  rename(MDD_Start = Start,
         MDD_End = End)

# check no duplicate Gene names
antidep_results$Gene %>%
  duplicated() %>%
  any()

# which genes are duplicated, look at entire row
antidep_results[duplicated(antidep_results$Gene) | duplicated(antidep_results$Gene, fromLast = TRUE), ] %>%
  arrange(Gene)

# Arrange so that genes with most methods and in MDD GWAS methods appear first
first_col <- length(c("ENSID", "Gene", "Chr", rep(c("Start", "End"),3))) + 1
last_col <- ncol(antidep_results) - 2
antidep_results <- antidep_results %>%
  arrange(desc(rowSums(.[first_col:last_col])))
write.csv(antidep_results, here::here("manuscript/tables/across_methods_and_mdd_gwas.csv"), row.names = F, quote = F)

# ------------------------------
# Save table of genes that are in antidep GWAS and may or may not be in MDD GWAS
write.csv(
  antidep_results %>%
    filter_at(vars(starts_with("mBAT_combo"), "SuSiEx", ends_with("SMR")), any_vars(. == TRUE)) %>%
    mutate(across(where(is.logical), ~ replace_na(., FALSE))), # change NA to FALSE, as that meant it was not identified in MDD GWAS as a high confidence gene
  here::here("manuscript/tables/across_methods_and_mdd_gwas_antidep_subset.csv"),
  row.names = F, quote = F
)

antidep_results %>%
  filter_at(vars(starts_with("mBAT_combo"), "SuSiEx", ends_with("SMR")), any_vars(. == TRUE)) %>%
  mutate(across(where(is.logical), ~ replace_na(., FALSE))) %>%
  pull(Gene) %>%
  duplicated(.) %>%
  sum()
# no duplicated genes

colname_descriptions <- c(
  "ENSID" = "Gene ID",
  "Gene" = "Gene name",
  "Chr" = "Chromosome",
  "mBATcombo_SuSiEx_Start" = "Start position of gene used in mBATcombo and SuSiEx",
  "mBATcombo_SuSiEx_End" = "End position of gene used in mBATcombo and SuSiEx",
  "SMR_Start" = "Start position of gene used in SMR",
  "SMR_End" = "End position of gene used in SMR",
  "Start_MDD" = "Start position of gene identified in MDD GWAS (Adams et al. 2025)",
  "End_MDD" = "End position identified in MDD GWAS (Adams et al. 2025)",
  "mBAT_combo_EUR_N06A" = "Gene identified with mBAT-combo in antidepressant GWAS meta-analysis (N06A) (EUR ancestry)",
  "mBAT_combo_EUR_N06AA" = "Gene identified with mBAT-combo in antidepressant GWAS meta-analysis (N06AA) (EUR ancestry)",
  "mBAT_combo_EUR_N06AB" = "Gene identified with mBAT-combo in antidepressant GWAS meta-analysis (N06AB) (EUR ancestry)",
  "SuSiEx" = "Gene identified in antidepressant GWAS meta-analysis (SuSiEx)",
  "blood_trait_eSMR" = "Gene identified with SMR in blood expression data in antidepressant GWAS meta-analysis (N06A) (EUR ancestry)",
  "blood_trait_mSMR" = "Gene identified with SMR in blood methylation data in antidepressant GWAS meta-analysis (N06A) (EUR ancestry)",
  "blood_trait_pSMR" = "Gene identified with SMR in blood proteomics data in antidepressant GWAS meta-analysis (N06A) (EUR ancestry)",
  "brainmeta_trait_eSMR" = "Gene identified with SMR in brain expression data in antidepressant GWAS meta-analysis (N06A) (EUR ancestry)",
  "brainmeta_trait_mSMR" = "Gene identified with SMR in brain methylation data in antidepressant GWAS meta-analysis (N06A) (EUR ancestry)",
  "brainmeta_trait_sSMR" = "Gene identified with SMR in brain splicing data in antidepressant GWAS meta-analysis (N06A) (EUR ancestry)",
  "MDD_GWAS_Nearest_gene" = "Gene identified with nearest gene method in MDD GWAS (Adams et al. 2025)",
  "MDD_GWAS_Fine_mapping" = "Gene identified with fine mapping method in MDD GWAS (Adams et al. 2025)",
  "MDD_GWAS_Expression" = "Gene identified with expression method in MDD GWAS (Adams et al. 2025)",
  "MDD_GWAS_Protein" = "Gene identified with protein method in MDD GWAS (Adams et al. 2025)",
  "MDD_GWAS_fastBAT" = "Gene identified with fastBAT method in MDD GWAS (Adams et al. 2025)",
  "MDD_GWAS_HMAGMA" = "Gene identified with HMAGMA method in MDD GWAS (Adams et al. 2025)",
  "MDD_GWAS_PsyOPS" = "Gene identified with PsyOPS method in MDD GWAS (Adams et al. 2025)"
)

# Write a .cols file
create_cols_meta(
  "manuscript/tables/across_methods_and_mdd_gwas_antidep_subset.csv",
  create_cols_meta,
  colname_descriptions
)
