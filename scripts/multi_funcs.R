# Helper functions used in multi.Rmd

# Plot ancestry PCs: eg. plot_ancestry_PCs("antidep-2501-mrmega-N06A.log")
plot_ancestry_PCs <- function(log_file_name){
  mrmega_log <- read_lines(here::here("meta", log_file_name))
  
  pcs_start_idx <- which(mrmega_log == "Principal components:") + 1
  pcs_end_idx <- which(mrmega_log == "Analysis finished.") - 1
  
  pcs <- read_table(str_c(mrmega_log[seq(pcs_start_idx, pcs_end_idx)], sep="\\n")) |>
    mutate(dataset = str_remove(PCs, fixed(".txt.gz"))) |>
    separate(dataset, sep = "-", into = c("cohort", "pheno", "cluster", "version"), remove = FALSE)
  
  plot <- ggplot(pcs, aes(x = PC0, y = PC1, label = str_glue("{cohort} [{cluster}]"))) +
    geom_point() +
    geom_text_repel()
  
  return(plot)
}

# Get sumstats. eg. get_sumstats("antidep-2501-mrmega-N06A.gz")
get_sumstats <- function(gz_file_name){
  mrmega <- read_tsv(here::here("meta", gz_file_name))
  return(mrmega)
}

# Split into a list. eg. get_assoc(mrmega, "N06A")
get_assoc <- function(mrmega, pheno_name){
  mrmega_assoc <- list(
     select(mrmega, ID = MarkerName, CHROM = Chromosome, POS = Position, P = `P-value_association`),
     select(mrmega, ID = MarkerName, CHROM = Chromosome, POS = Position, P = `P-value_ancestry_het`),
     select(mrmega, ID = MarkerName, CHROM = Chromosome, POS = Position, P = `P-value_residual_het`))
  
  names(mrmega_assoc) <- c(paste0(pheno_name, " Assoc"), 
                           paste0(pheno_name, " Ancestry"),
                           paste0(pheno_name, " Residual"))
  return(mrmega_assoc)
}

# Make manhattan plot, eg. plot_manhat(mrmega_assoc)
plot_manhat <- function(mrmega_assoc){
  plot <- manhattan(mrmega_assoc, legend_labels = names(mrmega_assoc), ntop = 1, sign_thresh = 5e-08, build = 38)
  return(plot)
}


# Get clumps. eg. get_clumps("antidep-2501-mrmega-N06A.clumps")
get_clumps <- function(clump_file_name){
  clumps <- read_tsv(here::here("meta", clump_file_name))
  
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


# Look up in GWAS catalogue eg. look_up_snps(clumps_gw_nomhc_gr, gwcat)
# gwcat <- get_cached_gwascat()
look_up_snps <- function(clumps_gw_nomhc_gr, gwcat){
  open_gwas <- phewas(variants = clumps_gw_nomhc_gr |> as_tibble() |> pull(SNP), pval=5e-8)
  
  gwcat_snps <-
    gwcat |>
    select(PUBMEDID, `DISEASE/TRAIT`, SNPS, MAPPED_TRAIT) |>
    separate_wider_delim(SNPS, delim = "; ", names_sep = "_", too_few = "align_start") |>
    pivot_longer(starts_with("SNP"), values_to = 'SNP') |>
    filter(!is.na(SNP)) |>
    select(-name)
  
  nomhc_snps <- clumps_gw_nomhc_gr |> as_tibble() |> pull(SNP)
  
  gwcat_snps |>
    filter(SNP %in% nomhc_snps) |>
    transmute(`DISEASE/TRAIT` = str_sub(`DISEASE/TRAIT`, 1, 50), SNP) |>
    arrange(SNP)
}
