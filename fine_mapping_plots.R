# Plot summary and locus zoom plots from SuSiEx output

#if (!require("susiexr", quietly = TRUE)) { devtools::install_github("ameliaes/susiexr") }
library(cowplot)
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)
library(stringr)
library(susiexR)

# ---- Read in variables:
plotsDir <- "fineMapping/plots/"

path_to_susiex_results <- "fineMapping/output"

ancestries <- fread('reference/ukb-ld-ref_ancestry.id') %>%
                pull(3) %>%
                unique()

bfile_paths <- path_to_susiex_results

sumstats_path <- path_to_susiex_results

sumstatsPathClean <- list.files(sumstats_path, pattern = ".ma$", full.names = TRUE)

sumstats <- lapply(sumstatsPathClean, function(sumstats_path){ 
  fread(sumstats_path)
  })

extract_ancestry_from_ma_path <- function(path) {
  str_remove(path, "\\.ma$") %>% # extract everything before the .ma file extension
    str_remove(., "^.*-") # extract the string after the last hyphen
}

ancestry_from_ma_path <- sapply(sumstatsPathClean, extract_ancestry_from_ma_path) %>%
  unname()

# check "ancestries" matches "ancestry_from_ma_path", else return error "Ancestries do not match"
if (!setequal(ancestries, ancestry_from_ma_path)) {
  stop("Ancestries are not defined in the file name of sumstats. 
       Please check documentation for information on file naming of sumstats.")
}

names(sumstats) <- ancestry_from_ma_path

# ---- Read in susiex results
results <- susiexR::format_results(path_to_susiex_results, ancestries)

# Check that each processed file type has the same number of fine mapped regions
nrow(results$summary) == length(results$cs) && length(results$cs) == length(results$snp)

# ---- Plot the relationship between:                                                       
# CS_LENGTH - number of SNPs in the credible set.                                            
# CS_PURITY - purity of the credible set.                                                    
# MAX_PIP - Maximum posterior inclusion probability (PIP) in the credible set..              
                                                                                          
png(paste0( plotsDir , "length_purity_maxPIP.png"), width = 1000, height = 600, res = 150)                     
print(plotPurityPIP(results$summary))                                                    
dev.off()                                                                                  
                                                                                          
# ---- Plot the probability the top SNP in the credible set is causal in each ancestry.                                                                                                  
png(paste0( plotsDir , "POST-HOC_PROB_POP.png"), width = 1000, height = 800, res = 150)                       
print(plotAncestryCausal(results$summary, ancestries = ancestries))                     
dev.off()                                                                                  
                                                                                          
# ---- Locus zoom plots
# Loop over fine mapped regions, creating a new plot file for each containing:
# Plot a region plot of -log10(p) vs SNP for each ancestry (1 ancestry plot per row)
# Another option is to overlay ancestries on top of each other and not colour by R2, but that might be messy?
# Bottom row is PIP vs SNP plot (the PIP is the overall PIP from susieX, a combination of all ancestries)s

# Define fine mapped regions:
# These will be each element of results$snp and results$cs
# They should be in the same order

# load in function for nextflow (todo: update function in susiexR package)
source( 'mainPlotNextFlow.R' )

for( i in 1:length(results$cs) ){
  main_plot <- mainPlotNextFlow(cs_results = results$cs[[i]], 
                      snp_results = results$snp[[i]],
                      sumstats = sumstats,
                      ancestries = ancestries,
                      plink2_path="plink2", 
                      bfile_paths = bfile_paths
                      )

  png(paste0(plotsDir, "region_plot_", main_plot$region, ".png"), width = 2300, height = 2300, res = 300)
  print(main_plot$plot)
  dev.off()
}

