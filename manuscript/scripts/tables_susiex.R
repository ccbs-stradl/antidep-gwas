# Create a table with significant results from the susiex analysis
# ---------------------
library(susiexR)
library(data.table)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(readr)

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
colname_descriptions <- c(
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


colname_descriptions_table <- tibble(column = names(colname_descriptions), description = colname_descriptions)

if(any(colname_descriptions_table$column != colnames(significant_snps_summary))){
  stop(glue("Column names in {file_name_summary} are not all described in colname_descriptions"))
}

write_tsv(colname_descriptions_table, here::here(glue('{file_name_summary}.cols')))

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
rm(colname_descriptions_table)
rm(colname_descriptions)

colname_descriptions <- c(
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

colname_descriptions_table <- tibble(column = names(colname_descriptions), description = colname_descriptions)

if(any(colname_descriptions_table$column != colnames(significant_snps_cs))){
  stop(glue("Column names in {file_name_cs} are not all described in colname_descriptions"))
}

write_tsv(colname_descriptions_table, here::here(glue('{file_name_cs}.cols')))

