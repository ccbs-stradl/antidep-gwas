# Create a summary table where each row has a method and the columns are the GWASs
# Value in cell is proportion of overlap with MDD GWAS
# ------------------------------
library(dplyr)

# Read in overlapping genes between antidep GWAS and MDD GWAS
#### NOTE MDD genes are from "supp_table_8B" which are high confidence genes. ####
antidep_results <- read.csv("manuscript/tables/across_methods_and_mdd_gwas_antidep_subset.csv")

# Create a col called MDD_GWAS_any
antidep_results <- antidep_results %>%
  mutate(MDD_GWAS_any = rowSums(across(starts_with("MDD_GWAS_")), na.rm = TRUE) > 0)


# Define a summary table with the columns: cross-ancestry, EUR_N06A, EUR_N06AA, EUR_N06AB
# and the rows: positional_mapping, fine_mapping, SMR
col_names <- c("cross-ancestry_N06A", "EUR_N06A", "EUR_N06AA", "EUR_N06AB")
row_names <- c("positional_mapping", "fine_mapping", "blood_trait_eSMR", 
               "blood_trait_mSMR", "blood_trait_pSMR", "brainmeta_trait_eSMR", 
               "brainmeta_trait_mSMR", "brainmeta_trait_sSMR")
summary_table <- data.frame(matrix(ncol = length(col_names), nrow = length(row_names)))
colnames(summary_table) <- col_names
rownames(summary_table) <- row_names


# Map row_names to columns in antidep_results
mapped_values <- list(
                  list("antidep_gwas" = "mBAT_combo_EUR_N06A",
                       "mdd_gwas" = "MDD_GWAS_fastBAT",
                       "summary_table_col" = "EUR_N06A",
                      "summary_table_row" = "positional_mapping"
                      ),
                  list("antidep_gwas" = "mBAT_combo_EUR_N06AA",
                       "mdd_gwas" = "MDD_GWAS_fastBAT",
                       "summary_table_col" = "EUR_N06AA",
                       "summary_table_row" = "positional_mapping"
                      ),
                  list("antidep_gwas" = "mBAT_combo_EUR_N06AB",
                       "mdd_gwas" = "MDD_GWAS_fastBAT",
                       "summary_table_col" = "EUR_N06AB",
                       "summary_table_row" = "positional_mapping"
                      ),
                  list(
                      "antidep_gwas" = "SuSiEx",
                      "mdd_gwas" = "MDD_GWAS_Fine_mapping",
                      "summary_table_col" = "cross-ancestry_N06A",
                      "summary_table_row" = "fine_mapping"
                      ),
                  list("antidep_gwas" = "blood_trait_eSMR",
                       "mdd_gwas" = "MDD_GWAS_any",
                       "summary_table_col" = "EUR_N06A",
                       "summary_table_row" = "blood_trait_eSMR"
                      ),
                  list("antidep_gwas" = "blood_trait_mSMR",
                       "mdd_gwas" = "MDD_GWAS_any",
                       "summary_table_col" = "EUR_N06A",
                       "summary_table_row" = "blood_trait_mSMR"
                      ),
                  list("antidep_gwas" = "blood_trait_pSMR",
                       "mdd_gwas" = "MDD_GWAS_any",
                       "summary_table_col" = "EUR_N06A",
                       "summary_table_row" = "blood_trait_pSMR"
                      ),
                  list("antidep_gwas" = "brainmeta_trait_eSMR",
                       "mdd_gwas" = "MDD_GWAS_any",
                       "summary_table_col" = "EUR_N06A",
                       "summary_table_row" = "brainmeta_trait_eSMR"
                      ),
                  list("antidep_gwas" = "brainmeta_trait_mSMR",
                       "mdd_gwas" = "MDD_GWAS_any",
                       "summary_table_col" = "EUR_N06A",
                       "summary_table_row" = "brainmeta_trait_mSMR"
                      ),
                  list("antidep_gwas" = "brainmeta_trait_sSMR",
                       "mdd_gwas" = "MDD_GWAS_any",
                       "summary_table_col" = "EUR_N06A",
                       "summary_table_row" = "brainmeta_trait_sSMR"
                      )
                  )

##################################
### FUNCTIONS ####################
##################################

# Define function to get the proportion and % of overlaping genes with MDD GWAS
get_cell_value <- function(antidep_col, mdd_col){
  antidep_genes_n <- antidep_results %>%
    filter(!!sym(antidep_col)) %>%
    nrow()

  mdd_genes_n <- antidep_results %>%
    filter(!!sym(mdd_col)) %>%
    nrow()

  antidep_and_mdd_genes_n <- antidep_results %>%
    filter(!!sym(antidep_col)) %>%
    filter(!!sym(mdd_col)) %>%
    nrow()

  total_genes <- 19000 # x-chromosome genes excluded

  # Perform pypher test to test whether the overlap between two gene sets
  # is more than what you’d expect by chance
  p_value <- phyper(antidep_and_mdd_genes_n - 1,
                    antidep_genes_n,
                    total_genes - antidep_genes_n,
                    mdd_genes_n,
                    lower.tail = F)

  # Round the p-value to 3 decimal places
  p_value_rounded = signif(p_value, digits = 3)

  # Get the proportion and percentage of overlapping genes
  proportion_gene_overlap <- paste0( antidep_and_mdd_genes_n , "/" , antidep_genes_n)
  percentange_gene_overlap <- paste0( round(antidep_and_mdd_genes_n / antidep_genes_n * 100, 2), "%")

  cell_value <- paste0(proportion_gene_overlap, " (", percentange_gene_overlap, "; p-value=", p_value_rounded, ")")

  return(cell_value)
}


# Convert the cell value to a dataframe with correct row name and col name
convert_to_dataframe <- function(antidep_col, mdd_col, row_name, col_name){
  cell_value <- get_cell_value(antidep_col, mdd_col)
  cell_value_df <- as.data.frame(cell_value)
  rownames(cell_value_df) <- row_name
  colnames(cell_value_df) <- paste0(colnames(cell_value_df), "_", col_name)
  return(cell_value_df)
}

main <- function(){

  # Fill in the summary table
  for(i in 1:length(mapped_values)){
    antidep_col <- mapped_values[[i]]$antidep_gwas
    mdd_col <- mapped_values[[i]]$mdd_gwas
    row_name <- mapped_values[[i]]$summary_table_row
    col_name <- mapped_values[[i]]$summary_table_col

    cell_value_df <- convert_to_dataframe(antidep_col, mdd_col, row_name, col_name)
    summary_table[row_name, col_name] <- cell_value_df
  }

  # Write the table to a csv file
  write.csv(summary_table, "manuscript/tables/antidep_gwas_mdd_gwas_summary_table.csv", row.names = T, quote = F)

  return(summary_table)

}

main()


# ------------
# Get rows that are in both MDD and antidep GWAS
# to get gene name from to insert into results section
# we could save these tables as extra sheets in Supp Table 8, as filtered views?

get_gene_rows <- function(antidep_col, mdd_col){

  antidep_results %>%
    filter(!!sym(antidep_col)) %>%
    filter(!!sym(mdd_col)) 

}

lapply(1:length(mapped_values), function(i) {
  antidep_col <- mapped_values[[i]]$antidep_gwas
  mdd_col <- mapped_values[[i]]$mdd_gwas

  
  get_gene_rows(antidep_col, mdd_col)
})
  




