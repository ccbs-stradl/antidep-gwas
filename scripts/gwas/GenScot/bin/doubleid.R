#!/usr/bin/env Rscript

# convert a FID/IID file to a double-id (IID/IID) file

library(dplyr)
library(readr)

args <- commandArgs(TRUE)

input_path <- args[1]
output_path <- args[2]

input <- read_table(input_path, col_types=cols(FID=col_character(), IID=col_character()))

output <- 
input |>
mutate(FID=IID) |>
select(FID, IID, everything())

write_tsv(output, output_path)
