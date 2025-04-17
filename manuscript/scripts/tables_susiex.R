# Create a table with significant results from the susiex analysis
# ---------------------
library(susiexR)
library(data.table)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)

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

write.csv(significant_snps_summary, here::here('manuscript/tables/susiex_significant_summary.csv'), row.names = FALSE, quote = FALSE)


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

write.csv(significant_snps_cs, here::here('manuscript/tables/susiex_significant_cs.csv'), row.names = FALSE, quote = FALSE)

# ---------------------
# .snp files
results$snp[items_to_keep] # this is a big file with all the SNPs, not sure if it's necessary to put in supplementary?

