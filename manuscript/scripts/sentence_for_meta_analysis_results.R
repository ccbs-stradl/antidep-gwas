# Text for meta-analysis result section
# "For the fixed effect meta-analyses of N06A the power equivalent to a
# case-control study for each ancestry was:
# AFRcases = X, AFRcontrols = ; EURcases = X, EURcontrols = …."

library("rjson")

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
  paste0(cluster, " cases = ",N , ", ", cluster, " controls = ", N, "; ")
}) %>% do.call(paste0, .)
