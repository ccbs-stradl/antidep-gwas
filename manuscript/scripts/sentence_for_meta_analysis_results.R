# Text for meta-analysis result section
# "For the fixed effect meta-analyses of N06A the power equivalent to a
# case-control study for each ancestry was:
# AFRcases = X, AFRcontrols = ; EURcases = X, EURcontrols = …."

library("rjson")
library(stringr)

json_files <- list.files("vcf/meta/GRCh38",
           pattern = "N06A-.*\\.json$",
           full.names = TRUE)

json_data <- lapply(json_files, function(json_file) {
  json <- fromJSON(file=json_file)
  as.data.frame(json)
})

lapply(json_data, function(json){
  cluster <- json[["cluster"]]
  neff <- json[["neff"]]
  N <- round(neff/2, 3)
  paste0(cluster, " = ",N , "; ")
}) %>% do.call(paste0, .)

# -----------------------------------------
# For the MR-MEGA numbers read in Neff from the cohort summary tables 
# (can also do this for fixed as it will return the same number)

summary_table_path <- here::here("manuscript", "tables", "meta_analysis_cohort_summary_table_mrmega.csv")
summary_table <- read.csv(summary_table_path)
summary_table_neff <- summary_table %>%
  rowwise() %>%
  mutate(total_neff = round(sum(c_across(ends_with("neff")), na.rm = TRUE)/2),3) %>%
  ungroup() %>%
  select(meta_analysis,cluster, total_neff) %>%
  filter(str_detect(meta_analysis,  "N06A$")) 

lapply(1:nrow(summary_table_neff), function(i){
  cluster <- summary_table_neff$cluster[i]
  N <- summary_table_neff$total_neff[i]
  paste0(cluster, " = ", N , "; ")
}) %>% do.call(paste0, .)

