# Create main table for mBat-combo results

library(data.table)
library(dplyr)
library(stringr)

# Change dir to "./mapsDir"
# for now use output from: "/Volumes/GenScotDepression/users/amelia/AMBER/maps_hg19"


# list.file("./mapsDir", full.names = TRUE, pattern = ".mbat.tsv")
mbat_output_paths <- list.files("/Volumes/GenScotDepression/users/amelia/AMBER/maps_hg19",
                          full.names = TRUE, pattern = ".mbat.tsv")

mbat_output <- lapply(mbat_output_paths, fread)

# Add ancestry to each data.table
names(mbat_output) <- mbat_output_paths %>%
                        str_extract("EUR|AFR|EAS|SAS|AMR") # change to ancestries used in meta-analysis

head(mbat_output)

# Remove items from list with 0 rows
keep_items <- sapply(mbat_output, function(x) nrow(x) > 0)
mbat_output <- mbat_output[keep_items]

# See which ancestries have any results:
names(mbat_output) # only EUR

# Filter table to Bonferroni corrected P values < 0.05
pcorrect_table <- function(mbat_table){
  mbat_table %>%
    mutate(P_mBAT_bonferroni = p.adjust(P_mBAT, method = "bonferroni")) %>%
    filter(P_mBAT_bonferroni < 0.05) %>%
    arrange(desc(Chisq_mBAT))
}

mbat_output_corrected <- lapply(mbat_output, pcorrect_table)

# Create main table
save_table <- function(i){
  ancestry <- names(mbat_output_corrected[i])
  write.csv(mbat_output_corrected[[i]], paste0("manuscript/tables/mBAT-combo.", ancestry, ".csv") , row.names = F, quote = F)
}

lapply(1:length(mbat_output_corrected), save_table)

# Get number of genes
lapply(mbat_output_corrected, nrow)

