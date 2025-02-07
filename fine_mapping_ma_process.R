library(data.table)
library(dplyr)

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# The first argument is the path to .ma
ma_path <- args[1]

# Read in sumsstats and remove rows where there are duplicate SNPs
sumstats <- fread( ma_path ) %>%
              distinct(SNP, .keep_all = TRUE) # keeps first distinct SNP row; .keep_all =T returns all columns

# Rewrite sumstats with suffix "noZero"
fwrite(sumstats, ma_path, sep = "\t")
