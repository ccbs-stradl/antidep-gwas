# Compare results with MDD GWAS (Adams et al. 2025)
# create a table with each row is a gene, and columns for different methods,
# with further columns for if the gene was identified in MDD GWAS (noting the method)
# ------------------------------
# Read in Supplementary Table 8 from Adams et al. 2025 (Gene mapping results summarized across methods, related to STAR Methods)
# https://www.cell.com/cell/fulltext/S0092-8674(24)01415-6

library(tidyr)
library(dplyr)
library(data.table)


mdd_genes <- fread("manuscript/MDD_GWAS_supp_table_8B.csv") %>%
  # add prefix "MDD_GWAS_" to columns: c(Nearest_gene Fine_mapping Expression Protein fastBAT HMAGMA PsyOPS)
  rename_with(~paste0("MDD_GWAS_", .), c(Nearest_gene, Fine_mapping, Expression, Protein, fastBAT, HMAGMA, PsyOPS)) %>%
  # rename Chrom_b37 pos_start_b37 pos_end_b37 to Chr     Start       End
  rename(Chr = Chrom_b37, Start = pos_start_b37, End = pos_end_b37) %>%
  # remove cols ending in b38
  select(-ends_with("b38")) %>%
  mutate(across(c(Chr, Start, End), ~as.character(.)))


# ------------------------------
# Read in mBAT-combo results
mBAT_combo_results <- fread("manuscript/tables/mBAT-combo.csv")

mBAT_combo_results %>%
  group_by(phenotype, ancestry) %>%
  summarise(n_genes = n_distinct(Gene),
            n_genes_in_MDD = sum(Gene %in% mdd_genes$ENSID))


# pivot from long table to wide table for phenotype
mBAT_combo_wide <- mBAT_combo_results %>%
  select(phenotype, Gene, gene_name, Chr, Start, End) %>%
  mutate(value = TRUE) %>%  # Create a column to hold TRUE values
  pivot_wider(
    names_from = phenotype,
    values_from = value,
    names_prefix = "mBAT_combo_EUR_",
    values_fill = list(value = FALSE)  # Fill missing values with FALSE
  )

# ------------------------------
# Read in susiex results
susiex_results <- data.frame(gene_name = c("PLCL2", "WNT3"),
                             SuSiEx = TRUE)

# Join with mbat combo results
antidep_results <- merge(mBAT_combo_wide, susiex_results, by = c("gene_name"), all = TRUE) %>%
  mutate(SuSiEx = ifelse(is.na(SuSiEx), FALSE, SuSiEx)) %>%
  mutate(across(starts_with("mBAT_combo"), ~ifelse(is.na(.), FALSE, .)))

# Check any NA
antidep_results %>%
  filter(is.na(Gene))

# Manually fill in gene details for WNT3 (https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000108379;r=17:44839872-44910520)
antidep_results[antidep_results$gene_name == "WNT3", c("Gene", "Chr", "Start", "End")] <- c("ENSG00000108379", "17", 44839872, 44910520 )

# ------------------------------
# Join with MDD GWAS results
antidep_results <- antidep_results %>%
  rename(ENSID = Gene, Gene = gene_name) %>%
  # convert Chr, Start and End cols to numeric values
  mutate(across(c(Chr, Start, End), ~as.character(.))) %>%
  merge(mdd_genes, by =c("ENSID", "Gene", "Chr"), all = TRUE)

# Some discrepancies between Start.x and Start.y and End.x and End.y
# Merge them into one col each keeping *.x if there is a difference in numeric value
# TODO

# Arrange so that genes with most methods and in MDD GWAS methods appear first
antidep_results <- antidep_results %>%
  arrange(desc(rowSums(.[6:16])))


write.csv(antidep_results, "manuscript/tables/across_methods_and_mdd_gwas.csv", row.names = F, quote = F)

# ------------------------------
# Save table of genes that are in antidep GWAS and may or may not be in MDD GWAS
write.csv(antidep_results %>%
            filter_at(vars(starts_with("mBAT_combo"), "SuSiEx"), any_vars(. == TRUE)) %>%
            rename(Start = Start.x, End = End.x) %>%
            select(-c(Start.y, End.y)) %>%
            mutate(across(everything(), ~ replace_na(., FALSE))), # change NA to FALSE, as that meant it was not identified in MDD GWAS as a high confidence gene
          "manuscript/tables/across_methods_and_mdd_gwas_antidep_subset.csv", row.names = F, quote = F)

