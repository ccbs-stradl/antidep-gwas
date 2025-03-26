get_clumps <- function(clump_file_name, mrmega){
  clumps <- read_tsv(here::here("meta", "antidep-2501", clump_file_name))
  
  names(clumps)[names(clumps) == "P-value_association"] <- "P"
  
  clumps_gw <- clumps |>
    filter(P <= 5e-8) |>
    left_join(mrmega, by = c(`#CHROM` = "Chromosome", "POS" = "Position"))
  
  clumps_gw_gr <- clumps_gw |>
    select(seqnames = `#CHROM`, start = POS, width = 1, P, SNP = MarkerName) |>
    as_granges() |>
    set_genome_info(genome = 'hg38')
  
  mhc_gr <- tibble(seqnames = 6, start = 28510120, end = 33480577) |>
    as_granges() |>
    set_genome_info(genome = 'hg38')
  
  clumps_gw_mhc_gr <-
    clumps_gw_gr |>
    filter_by_overlaps(mhc_gr)
  
  clumps_gw_nomhc_gr <-
    clumps_gw_gr |>
    filter_by_non_overlaps(mhc_gr)
  
  return(clumps_gw_nomhc_gr)
}

clumps_gw_nomhc_gr_N06A <- get_clumps("antidep-2501-mrmega-N06A-DIV.clumps.tsv", mrmega_n06a_assoc)