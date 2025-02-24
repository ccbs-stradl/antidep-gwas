# Compare results with MDD GWAS (Adams et al. 2025)
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

genes_at_least_one_method <- fread("manuscript/MDD_GWAS_supp_table_8A.csv")
genes_high_conf <- fread("manuscript/MDD_GWAS_supp_table_8B.csv")

# ------------------------------
# Read in mBAT-combo results
mBAT_combo_results <- fread("manuscript/tables/mBAT-combo.csv")

# Identify which genes overlap with MDD GWAS and from which methods
mBAT_combo_results %>%
  filter(mBAT_combo_results$Gene %in% genes_at_least_one_method$ENSID)

# idea for figure? (all EUR)
mBAT_combo_results %>%
  filter(mBAT_combo_results$Gene %in% genes_high_conf$ENSID) %>%
  arrange(P_mBATcombo) %>%
  ggplot(data = .) +
  geom_point(aes(y = gene_name,
                 x = Chisq_mBAT,
                 colour = phenotype))

# For each ancestry and antidep subgroup combo are any of those genes in MDD GWAS?
mBAT_combo_results %>%
  group_by(phenotype, ancestry) %>%
  summarise(n_genes = n_distinct(Gene),
            n_genes_in_MDD = sum(Gene %in% genes_at_least_one_method$ENSID))

mBAT_combo_results %>%
  group_by(phenotype, ancestry) %>%
  summarise(n_genes = n_distinct(Gene),
            n_genes_in_MDD = sum(Gene %in% genes_high_conf$ENSID))



# ------------------------------
# SuSiEx

# ------------------------------
