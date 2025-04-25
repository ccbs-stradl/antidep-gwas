# Create main table for mBat-combo results

library(data.table)
library(dplyr)
library(stringr)
library(readr)

# Read in results from mBAT-combo
mbat_output_paths <- list.files(here::here("results", "maps", "mbat", "hg19"),
                          full.names = TRUE, pattern = ".mbat.tsv")

mbat_output <- lapply(mbat_output_paths, fread)

# Add ancestry and phenotype names to each data.table
names(mbat_output) <- basename(mbat_output_paths)

head(mbat_output)

# Remove items from list with 0 rows
keep_items <- sapply(mbat_output, function(x) nrow(x) > 0)
mbat_output <- mbat_output[keep_items]

# See which ancestries have any results:
names(mbat_output)

# Filter table to Bonferroni corrected P values < 0.05
pcorrect_table <- function(mbat_table){
  mbat_table %>%
    mutate(P_mBATcombo_bonferroni = p.adjust(P_mBATcombo, method = "bonferroni")) %>%
    filter(P_mBATcombo_bonferroni < 0.05) %>%
    arrange(desc(Chisq_mBAT))
}

mbat_output_corrected <- lapply(mbat_output, pcorrect_table)

# Add a new column of the phenotype and ancestry
for(i in 1:length(mbat_output_corrected)){
  split_name <- str_split(names(mbat_output_corrected)[i], "-")
  mbat_output_corrected[[i]]$phenotype <- split_name[[1]][4]
  mbat_output_corrected[[i]]$ancestry <- split_name[[1]][5]
  mbat_output_corrected[[i]] <- relocate(mbat_output_corrected[[i]], phenotype, ancestry)
}


# Get number of genes per ancestry/pheno
lapply(mbat_output_corrected, nrow)

# Combine into one table to save
mbat_results <- do.call(rbind, mbat_output_corrected)
file_name <- "manuscript/tables/mBAT-combo.csv"
write.csv(mbat_results, file_name , row.names = F, quote = F)

# ------------------------------------------
# Create the metadata sidecar .cols file 
# containing descriptions of column names


colname_descriptions <- c("phenotype" = "anti-depressant ATC code",
                          "ancestry" = "genetic ancestry",
                          "Gene" = "Ensembl gene ID",
                          "gene_name" = "gene name",
                          "Chr" = "chromosome",
                          "Start" = "start base position of the gene (Gr37)",
                          "End" = "end base position of the gene (Gr37)",
                          "No.SNPs" = "Number of SNPs in gene",
                          "SNP_start" = "rsID of SNP at atart of the gene",
                          "SNP_end" = "rsID of SNP at end of the gene",
                          "TopSNP" = "rsID of the top SNP in the gene",
                          "TopSNP_Pvalue" = "p-value of the top SNP in the gene",
                          "No.Eigenvalues" = "number of eigenvalues corresponding to the gamma parameter in mBAT",
                          "Chisq_mBAT" = "Chi-squared statistic for mBAT",
                          "P_mBATcombo" = "p-value for mBAT-combo",
                          "P_mBAT" = "p-value for mBAT",
                          "Chisq_fastBAT" = "Chi-squared statistic for fastBAT",
                          "P_fastBAT" = "p-value for fastBAT",
                          "P_mBATcombo_bonferroni" = "Bonferroni corrected p-value for mBAT-combo")

colname_descriptions_table <- tibble(column = names(colname_descriptions), description = colname_descriptions)

if(any(colname_descriptions_table$column != colnames(mbat_results))){
  stop(glue("Column names in {file_name} are not all described in colname_descriptions"))
}

write_tsv(colname_descriptions_table, here::here(glue('{file_name}.cols')))

# ------------------------------------------
# Create a subtable to show which genes are overlapping across phenotypes and ancestries
all_genes <- mbat_output_corrected %>%
  do.call(rbind,.) %>%
  pull(gene_name) %>%
  unique()

genes_df <- data.frame(gene_name = all_genes, matrix(ncol = length(names(mbat_output_corrected)), nrow = length(all_genes)))

split_name <- sapply(1:length(mbat_output_corrected), function(i) {
  split_name <- str_split(names(mbat_output_corrected[i]), "-")[[1]]
  paste0(split_name[4],"-", split_name[5])
})

colnames(genes_df)[-1] <- split_name

for(i in 1:length(mbat_output_corrected) ){
  genes_df[,1+i] <- genes_df$gene_name %in% mbat_output_corrected[[i]]$gene_name
}

overlapping <- genes_df %>%
  rowwise() %>%
  filter(sum(c_across(-gene_name)) >= 2) %>%
  ungroup()

overlapping

write.csv(overlapping, paste0("manuscript/tables/mBAT-combo_overlapping_genes.csv") , row.names = F, quote = F)




