# Create a phenogram diagram for overlap of MDD and antidep GWAS genes:
# https://visualization.ritchielab.org/phenograms/document



library(dplyr)
library(data.table)

# Read in overlapping genes between antidep GWAS and MDD GWAS
#### NOTE MDD genes are from "supp_table_8B" which are high confidence genes. ####
antidep_results <- read.csv("manuscript/tables/across_methods_and_mdd_gwas.csv")

head(antidep_results)

# system("wget https://visualization.ritchielab.org/downloads/phenogram-sample.txt")
example <- fread("phenogram-sample.txt")
system("head phenogram-sample.txt")
example

# Create a table with:
# (these should actually be lowercase)
# gene_name, Gene
# CHR
# POS	X	Base pair location of the SNP or in the case of a region, the starting location
# PHENOTYPE	X	Phenotype name (required unless plotting chromosomes only)
# END		Ending base pair location for a region
# NOTE or ANNOTATION		Locations on the plot will be annotated to the right of the chromosomes with values from this column. The labels can be no more than 10 characters in length.
# ETHNICITY or RACE or GROUP		Specifies the race/ethnicity for the result. When more than one occurs in the input, the plot will use different shapes to represent them.

# Convert table from wide to long
antidep_results %>%
  # if any of the columns starting with "MDD_GWAS" is TRUE, then the gene is in MDD GWAS
  # and new column MDD is TRUE, else FALSE
  mutate(MDD = rowSums(select(., starts_with("MDD_GWAS"))) > 0) %>%
  # same for antidep gwas, select starts with mBAT_combo and also SuSiEx column
  mutate(antidep = rowSums(select(., starts_with("mBAT_combo"), SuSiEx)) > 0) %>%
  # select(gene_name, Gene, CHR = Chr, POS = Start, END = End, MDD, antidep ) %>%
  # change from wide to long to create a phenotype column for "MDD" and "antidep" cols
  pivot_longer(cols = c("MDD", "antidep"), names_to = "phenotype", values_to = "GENE_IN_PHENOTYPE") %>%
  filter(GENE_IN_PHENOTYPE) %>%
  # Create POS and END cols, taking POS = Start.x, END = End.x if PHENOTYPE = antidep, and *.y if phenotype MDD
  mutate(pos = ifelse(phenotype == "antidep", Start.x, Start.y),
         end = ifelse(phenotype == "antidep", End.x, End.y)) %>%
  select(chrom = Chr, pos, phenotype) %>%
  # have to drop some MDD rows with NA in for now, todo fix!
  drop_na() %>%
  write.table(.,"manuscript/tables/phenogram_input.txt", row.names = F, quote = F, sep = "\t")



