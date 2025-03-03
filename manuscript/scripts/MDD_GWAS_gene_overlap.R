# Compare results with MDD GWAS (Adams et al. 2025)
# create a table with each row is a gene, and columns for different methods,
# with further columns for if the gene was identified in MDD GWAS (noting the method)
# ------------------------------
# Read in Supplementary Table 8 from Adams et al. 2025 (Gene mapping results summarized across methods, related to STAR Methods)
# https://www.cell.com/cell/fulltext/S0092-8674(24)01415-6
# Trans-ancestry genome-wide study of depression identifies 697 associations implicating cell types and pharmacotherapies
# Major Depressive Disorder Working Group of the Psychiatric Genomics Consortium1
# Comparison of gene mapping method. Table S8A lists genes identified by at least one method.
# There is a column for each method (Nearest_gene to PsyOPS) with TRUE if the gene was identified
# or prioritised by a given method (FALSE otherwise). Chromosome and gene start and end
# positions are listed for genome builds 37 (GRCh37/hg19) and 38 (GRCh38/hg38).
# High-confidence genes identified from finemapping, expression, or protein are also reported in the second sheet (B).

library(tidyr)

mdd_genes <- fread("manuscript/MDD_GWAS_supp_table_8A.csv")

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
antidep_results %>%
  mutate(MDD_GWAS = ifelse(gene_name %in% mdd_genes$Gene, TRUE, FALSE))

mdd_methods <- colnames(mdd_genes)[3:9]

add_mdd_cols <- function(mdd_method){
mdd_genes_tmp <- mdd_genes %>%
  filter(!!sym(mdd_methods[i]))

new_col_name <- paste0("MDD_GWAS_", mdd_methods[i])

antidep_results <- antidep_results %>%
  mutate(!!new_col_name := ifelse(gene_name %in% mdd_genes_tmp$Gene, TRUE, FALSE))

return(antidep_results)
}

for(i in 1:length(mdd_methods)){
  antidep_results <- add_mdd_cols(mdd_methods[i])
}

head(antidep_results)

# Arrange so that genes with most methods and in MDD GWAS methods appear first
antidep_results <- antidep_results %>%
  arrange(desc(rowSums(.[-(1:5)])))


write.csv(antidep_results, "manuscript/tables/across_methods_and_mdd_gwas.csv", row.names = F, quote = F)


# -------------------------------------
# Compare genes in antidep-gwas that ARE and ARE NOT in the MDD GWAS
antidep_results %>%
  filter_at(vars(starts_with("MDD_GWAS")), all_vars(. == FALSE)) %>%
  nrow()
# 97

antidep_results %>%
  filter_at(vars(starts_with("MDD_GWAS")), any_vars(. == TRUE)) %>%
  nrow()
# 170

head(antidep_results)
# ------------------------------
write.csv(filter_at(antidep_results, vars(starts_with("MDD_GWAS")), all_vars(. == FALSE)) %>% relocate(Gene),
          "manuscript/tables/across_methods_and_mdd_gwas_not_in_mdd_gwas.csv", row.names = F, quote = F)


write.csv(filter_at(antidep_results, vars(starts_with("MDD_GWAS")), any_vars(. == TRUE)) %>% relocate(Gene) ,
          "manuscript/tables/across_methods_and_mdd_gwas_in_mdd_gwas.csv", row.names = F, quote = F)

