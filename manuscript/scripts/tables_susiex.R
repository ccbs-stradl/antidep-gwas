# Create a table with significant results from the susiex analysis
# ---------------------
library(susiexR)
library(data.table)
library(dplyr)
library(glue)
library(tidyr)
library(purrr)
library(stringr)
library(readr)
library(openxlsx)
source(here::here("manuscript/scripts/supplementary_tables_excell_create_cols_meta_FUN.R"))
# source(here::here("manuscript/scripts/supplementary_tables_excel_functions.R"))

# ---------------------
# Read in the results from the susiex analysis
path_to_susiex_results <- here::here("results/fineMapping/output")

# Get the ancestries IN THE ORDER THEY WERE SUPPLIED TO SUSIEX
# this can be obtained from file name
ancestries_susiex_order <- list.files(path_to_susiex_results, pattern = ".summary")[1] %>% # the first one will do
  str_extract("(?<=\\.)([A-Za-z-]+)(?=\\.)") %>% # extract the ancestry from the file name
  str_split("-") %>% unlist()  # split the ancestry by hyphen


# format the results into a list of tables
results <- susiexR::format_results(path_to_susiex_results, ancestries_susiex_order)

# ---------------------
# .summary files
# Create a table with the significant results

significant_snps_summary <- results$summary %>%
                      filter(MAX_PIP > 0.8) %>%
                      relocate(CHR, BP_START, BP_END, .before = CS_ID) %>%
                      separate(REF_ALLELE, into = paste0("REF_ALLELE_", ancestries_susiex_order), sep = ",") %>%
                      separate(ALT_ALLELE, into = paste0("ALT_ALLELE_", ancestries_susiex_order), sep = ",") %>%
                      separate(REF_FRQ, into = paste0("REF_FRQ_", ancestries_susiex_order), sep = ",") %>%
                      separate(BETA, into = paste0("BETA_", ancestries_susiex_order), sep = ",") %>%
                      separate(SE, into = paste0("SE_", ancestries_susiex_order), sep = ",") %>%
                      separate(`-LOG10P`, into = paste0("-LOG10P_", ancestries_susiex_order), sep = ",") %>%
                      rename_with(~ reduce(seq_along(ancestries_susiex_order),
                                           .init = .x,
                                           ~ str_replace(.x, as.character(.y), paste0("_", ancestries_susiex_order[.y]))),
                                  starts_with("POST-HOC_PROB_POP"))
file_name_summary <- 'manuscript/tables/susiex_significant_summary.csv'
write.csv(significant_snps_summary, here::here(file_name_summary), row.names = FALSE, quote = FALSE)

# Write out column descriptions
colname_descriptions_summary <- c(
  "CHR" = "Chromosome",
  "BP_START" = "Start position of the finemapping region",
  "BP_END" = "End position of the finemapping region",
  "CS_ID" = "Credible Set ID",
  "CS_LENGTH" = "Size (number of SNPs) of the credible set",
  "CS_PURITY" = "Purity of the credible set",
  "MAX_PIP_SNP" = "SNP in the credible set that had the largest posterior inclusion probability (PIP)",
  "BP" = "The base pair coordinate of the MAX_PIP_SNP in Gr37",
  setNames(
    paste0("Reference allele of the MAX_PIP_SNP for ancestry: ", ancestries_susiex_order),
    paste0("REF_ALLELE_", ancestries_susiex_order)
  ),
  setNames(
    paste0("Alternate allele of the MAX_PIP_SNP for ancestry: ", ancestries_susiex_order),
    paste0("ALT_ALLELE_", ancestries_susiex_order)
  ),
  setNames(
    paste0("Reference allele frequency for ancestry: ", ancestries_susiex_order),
    paste0("REF_FRQ_", ancestries_susiex_order)
  ),
  setNames(
    paste0("Marginal per-allele effect size of the MAX_PIP_SNP with respect to the reference allele for ancestry: ", ancestries_susiex_order),
    paste0("BETA_", ancestries_susiex_order)
  ),
  setNames(
    paste0("Standard error of the marginal per-allele effect size of the MAX_PIP_SNP for ancestry: ", ancestries_susiex_order),
    paste0("SE_", ancestries_susiex_order)
  ),
  setNames(
    paste0("-log10 of the marginal p-value of the MAX_PIP_SNP for ancestry: ", ancestries_susiex_order),
    paste0("-LOG10P_", ancestries_susiex_order)
  ),
  "MAX_PIP" = "Maximum posterior inclusion probability (PIP) in the credible set.",
  setNames(
    paste0("Post hoc probability credible set manifest causal in population: ", ancestries_susiex_order),
    paste0("POST-HOC_PROB_POP_", ancestries_susiex_order)
  )
)


colname_descriptions_tables_summary <- tibble(column = names(colname_descriptions_summary), description = colname_descriptions_summary)

if(any(colname_descriptions_tables_summary$column != colnames(significant_snps_summary))){
  stop(glue("Column names in {file_name_summary} are not all described in colname_descriptions_summary"))
}

write_tsv(colname_descriptions_tables_summary, here::here(glue('{file_name_summary}.cols')))

items_to_keep <- sapply(results$cs, function(results_table){
  sig <- results_table %>%
    filter(CHR %in% significant_snps_summary$CHR & BP_START %in% significant_snps_summary$BP_START)
  nrow(sig) > 0
})

# ---------------------
# .cs files
significant_snps_cs <- results$cs[items_to_keep] %>%
  bind_rows() %>%
  relocate(CHR, BP_START, BP_END, .before = CS_ID) %>%
  separate(REF_ALLELE, into = paste0("REF_ALLELE_", ancestries_susiex_order), sep = ",") %>%
  separate(ALT_ALLELE, into = paste0("ALT_ALLELE_", ancestries_susiex_order), sep = ",") %>%
  separate(REF_FRQ, into = paste0("REF_FRQ_", ancestries_susiex_order), sep = ",") %>%
  separate(BETA, into = paste0("BETA_", ancestries_susiex_order), sep = ",") %>%
  separate(SE, into = paste0("SE_", ancestries_susiex_order), sep = ",") %>%
  separate(`-LOG10P`, into = paste0("-LOG10P_", ancestries_susiex_order), sep = ",")

file_name_cs <- 'manuscript/tables/susiex_significant_cs.csv'
write.csv(significant_snps_cs, here::here(file_name_cs), row.names = FALSE, quote = FALSE)

# Write out column descriptions
# clear environment

colname_descriptions_cs <- c(
  "CHR" = "Chromosome",
  "BP_START" = "Start position of the finemapping region",
  "BP_END" = "End position of the finemapping region",
  "CS_ID" = "Credible Set ID",
  "SNP" = "SNP identifier",
  "BP" = "The base pair coordinate of the SNP in Gr37",
  setNames(
    paste0("Reference allele of the SNP for ancestry: ", ancestries_susiex_order),
    paste0("REF_ALLELE_", ancestries_susiex_order)
  ),
  setNames(
    paste0("Alternate allele of the SNP for ancestry: ", ancestries_susiex_order),
    paste0("ALT_ALLELE_", ancestries_susiex_order)
  ),
  setNames(
    paste0("Reference allele frequency for ancestry: ", ancestries_susiex_order),
    paste0("REF_FRQ_", ancestries_susiex_order)
  ),
  setNames(
    paste0("Marginal per-allele effect size of the SNP with respect to the reference allele for ancestry: ", ancestries_susiex_order),
    paste0("BETA_", ancestries_susiex_order)
  ),
  setNames(
    paste0("Standard error for ancestry: ", ancestries_susiex_order),
    paste0("SE_", ancestries_susiex_order)
  ),
  setNames(
    paste0("-log10(p-value) for ancestry: ", ancestries_susiex_order),
    paste0("-LOG10P_", ancestries_susiex_order)
  ),
  "CS_PIP" = "Posterior inclusion probability (PIP) of the SNP in the credible set.",
  "OVRL_PIP" = "Posterior inclusion probability (PIP) of the SNP in the entire region."
)

colname_descriptions_table_cs <- tibble(column = names(colname_descriptions_cs), description = colname_descriptions_cs)

if(any(colname_descriptions_table_cs$column != colnames(significant_snps_cs))){
  stop(glue("Column names in {file_name_cs} are not all described in colname_descriptions_cs"))
}

write_tsv(colname_descriptions_table_cs, here::here(glue('{file_name_cs}.cols')))

# ---------------------
# Create an xlsx with all the non-null SuSiEx results to upload to FigShare
# one sheet for: summary, cs and snp

# Same formatting as above for summary, 
# but with PIP filter removed to get significant results
summary <- results$summary %>%
  relocate(CHR, BP_START, BP_END, .before = CS_ID) %>%
  separate(REF_ALLELE, into = paste0("REF_ALLELE_", ancestries_susiex_order), sep = ",") %>%
  separate(ALT_ALLELE, into = paste0("ALT_ALLELE_", ancestries_susiex_order), sep = ",") %>%
  separate(REF_FRQ, into = paste0("REF_FRQ_", ancestries_susiex_order), sep = ",") %>%
  separate(BETA, into = paste0("BETA_", ancestries_susiex_order), sep = ",") %>%
  separate(SE, into = paste0("SE_", ancestries_susiex_order), sep = ",") %>%
  separate(`-LOG10P`, into = paste0("-LOG10P_", ancestries_susiex_order), sep = ",") %>%
  rename_with(~ reduce(seq_along(ancestries_susiex_order),
                       .init = .x,
                       ~ str_replace(.x, as.character(.y), paste0("_", ancestries_susiex_order[.y]))),
              starts_with("POST-HOC_PROB_POP"))

# Same formatting as above for cs
# but without subsetting by [items to keep]
cs <- results$cs %>%
  bind_rows() %>%
  relocate(CHR, BP_START, BP_END, .before = CS_ID) %>%
  separate(REF_ALLELE, into = paste0("REF_ALLELE_", ancestries_susiex_order), sep = ",") %>%
  separate(ALT_ALLELE, into = paste0("ALT_ALLELE_", ancestries_susiex_order), sep = ",") %>%
  separate(REF_FRQ, into = paste0("REF_FRQ_", ancestries_susiex_order), sep = ",") %>%
  separate(BETA, into = paste0("BETA_", ancestries_susiex_order), sep = ",") %>%
  separate(SE, into = paste0("SE_", ancestries_susiex_order), sep = ",") %>%
  separate(`-LOG10P`, into = paste0("-LOG10P_", ancestries_susiex_order), sep = ",")

# results$snp
# replace "Pop1", "Pop2", "Pop3" etc with ancestry names from ancestries_susiex_order
# Build named vector: names are what to replace (Pop1, Pop2...), values are the ancestry names
pop_map <- setNames(ancestries_susiex_order, paste0("Pop", seq_along(ancestries_susiex_order)))

snp <- results$snp %>%
  bind_rows() %>%
  relocate(CHR, BP_START, BP_END, .before = BP) %>%
  rename_with(~ str_replace_all(.x, pop_map))

colname_descriptions_snp <- c(
  "CHR" = "Chromosome",
  "BP_START" = "Start position of the finemapping region",
  "BP_END" = "End position of the finemapping region",
  "BP" = "The base pair coordinate of the SNP in Gr37",
  "PIP(CS1)" = "The posterior inclusion probability (PIP) of the SNP in credible set 1.",
  "LogBF(CS1,AFR)" = "The Logarithm of Bayes factor for credible set 1, in AFR.",
  "LogBF(CS1,SAS)" = "The Logarithm of Bayes factor for credible set 1, in SAS.",
  "LogBF(CS1,EUR)" = "The Logarithm of Bayes factor for credible set 1, in EUR."
)

# Save csvs and their sidecar .cols
file_name_summary_all <- 'manuscript/tables/susiex_all_cs.csv'
file_name_cs_all <- 'manuscript/tables/susiex_all_summary.csv'
file_name_snp_all <- 'manuscript/tables/susiex_all_snp.csv'

write.csv(summary, here::here(file_name_summary_all), quote = TRUE, row.names = FALSE)
write.csv(cs, here::here(file_name_cs_all), quote = TRUE, row.names = FALSE)
write.csv(snp, here::here(file_name_snp_all), quote = TRUE, row.names = FALSE)

# Write out column descriptions
mapply(create_cols_meta, 
       file_name = c(file_name_summary_all, file_name_cs_all, file_name_snp_all),
       table_variable_name = c("summary", "cs", "snp"),
       colname_descriptions = list(colname_descriptions_summary, colname_descriptions_cs, colname_descriptions_snp)
)

# Create xlsx
table_index <- 0
create_table(
  paths = rep("manuscript/tables/", 3),
  regex = c("susiex_all_summary.csv",
            "susiex_all_cs.csv",
            "susiex_all_snp.csv"),
  sheet_names = c("SuSiEx summary",
                  "SuSiEx credible sets",
                  "SuSiEx SNPs"),
  excel_file_name = here::here(glue("manuscript/tables/online_results_SuSiEx.xlsx")),
  table_index,
  legend_title = "Non-null output from SuSiEx",
  legend_text_prefix = "Results are split into ",
  legend_text_sections = c("summary sheet with credible set level information",
                           "credible set sheet containing information for all the SNPs included in credible sets",
                           "snp sheet containing information for all the SNPs that are used in the fine-mapping algorithm. These are the SNPs that are located in the specified fine-mapping region, available in the GWAS summary statistics and the reference panel for at least one population, and survived the minor allele frequency filtering."),
  cell_title_width = 39,
  cell_title_height = 49
)
