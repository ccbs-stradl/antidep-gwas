# Copy main results from datastore to here so they can be git tracked.
# Create .cols sidecar meta data for SMR results
# Paste sentences to go into manuscript

# load libraries
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
library(readr)
library(glue)
library(here)

# Load function to create .cols sidecar meta data file
source(here::here("manuscript/scripts/supplementary_tables_excell_create_cols_meta_FUN.R"))

# -----------------------------------------------
# Copy results to this projects directory, from datastore
# First I made a sym link to the datastore dir (see docs/ for how to do this)
# Read in all SMR results for blood and brain
paths <- list(
  "results/maps/smr/blood/trait_eSMR.merged.tsv",
  "results/maps/smr/blood/trait_mSMR.merged.tsv",
  "results/maps/smr/blood/trait_pSMR.merged.tsv",
  "results/maps/smr/brainmeta/trait_eSMR.merged.tsv",
  "results/maps/smr/brainmeta/trait_mSMR.merged.tsv",
  "results/maps/smr/brainmeta/trait_sSMR.merged.tsv"
)

# Create col name descriptions, same for all SMR results
colname_descriptions <-
  c("Gene" = "gene name",
    "qtl_name" = "QTL name",
    "probeID" = "probe ID",
    "ProbeChr" = "probe chromosome",
    "Probe_bp" = "probe position",
    "topSNP" = "SNP name",
    "topSNP_chr" = "SNP chromosome",
    "topSNP_bp" = "SNP position",
    "A1" = "the effect (coded) allele",
    "A2" = "the other allele",
    "Freq" = "frequency of the effect allele (estimated from the reference samples)",
    "b_GWAS" = "effect size from GWAS",
    "se_GWAS" = "standard error from GWAS",
    "p_GWAS" = "p-value from GWAS",
    "b_eQTL" = "effect size from eQTL study",
    "se_eQTL" = "standard error from eQTL study",
    "p_eQTL" = "p-value from eQTL study",
    "b_SMR" = "effect size from SMR",
    "se_SMR" = "standard error from SMR",
    "p_SMR" = "p-value from SMR",
    "p_SMR_Bonferroni" = "Bonferroni corrected p-value from SMR",
    "p_HEIDI" = "p-value from HEIDI (HEterogeneity In Depedent Instruments) test",
    "nsnp_HEIDI" = "number of SNPs used in the HEIDI test",
    "gene_id" = "Gene ID",
    "chr" = "Chromosome",
    "start" = "Start position",
    "end" = "end position",
    "strand" = "strand",
    "GWAS_LOCUS" = "gwas locus",
    "Lead_SNP" = "lead SNP",
    "Lead_SNP_BP" = "lead SNP position"
  )

# -----------------------------------------------
# rename file from "blood/trait_eSMR.merged.tsv" to 
# "blood_trait_eSMR.merged.tsv" etc.
# then move this file to manuscripts/tables/
rename_file_and_move <- function(old_path, new_path_prefix){
  # rename file
  tissue_type <- basename(dirname(old_path)) # extracts the last directory name
  new_file_name <- str_c(tissue_type, '_', basename(old_path))
  
  # full new path
  new_path <- paste0(new_path_prefix, "/", new_file_name)
  
  # copy results and rename with new name
  file.copy(
    from = here::here(old_path),
    to = here::here(new_path),
    overwrite = TRUE
  )
  
  return(new_path)
}

# Read in results
read_results <- function(path){
  
  results <- fread(here::here(path))

  return(results)
}

# Extract names from path
extract_names <- function(path){
  # change "blood_trait_eSMR.merged.tsv" to "blood_eSMR"

  new_name <- basename(path) %>%
    str_remove_all("\\.merged\\.tsv") %>%
    str_replace_all("trait_", "") %>%
    str_replace("brainmeta", "brain") # Simplify "brainmeta" to "brain"

  return(new_name)
}

# Reorder rows in results
reorder_results <- function(tables_list){
  # Check if required columns exist
  sapply(tables_list, function(df) {
    if (!all(c("p_SMR", "p_HEIDI") %in% colnames(df))) {
      stop("Required columns 'p_SMR' and 'p_HEIDI' not found in the data frame.")
    }
  })
  
  # Reorder the rows in each table by the "p" column
  for (i in seq_along(tables_list)) {
    tables_list[[i]] <- tables_list[[i]] %>%
      arrange(desc(p_HEIDI)) %>%
      arrange(p_SMR) %>%
      mutate(p_SMR_Bonferroni = p.adjust(p_SMR, method = "bonferroni")) %>%
      relocate(
        p_SMR_Bonferroni,
        .after = p_SMR
      )
  }
  
  return(tables_list)
}


# Paste to console sentence to include in manuscript
paste_sentences <- function(tables_list){
  # Get significant results for each tissue and omic type
  gene_names <- lapply(1:length(tables_list), function(i) {
    genes <- tables_list[[i]] %>%
      filter(p_SMR_Bonferroni < 0.05 & p_HEIDI > 0.05)
    
    if (nrow(genes) > 0) {
      gene_names <- genes %>%
        pull(if_else("Gene" %in% colnames(genes), "Gene", "index"))
    } else {
      gene_names <- NULL
    }
    
    # Number of genes
    n_genes <- length(gene_names)
    
    # Name of the tissue and omic type
    tissue_omic <- names(tables_list)[i] %>%
      str_replace("_", " and ")
    
    # paste sentence
    message("", n_genes, " genes in ", tissue_omic)
    
    return(gene_names)
  })
  
  # Paste sentence on how many genes unique across all tissue and omic types
  n_unique_genes <- length(unique(unlist(gene_names)))
  message("Identified ", n_unique_genes, " unique genes across all tissue and omic types.")
  
}

main <- function(paths){
  new_paths <- lapply(paths, rename_file_and_move, new_path_prefix = "manuscript/tables")
  
  # Read in results
  tables_list <- lapply(new_paths, read_results)
  
  # Name results by their short name, eg. brain_eSMR
  names(tables_list) <- lapply(new_paths, extract_names)
  
  # Reorder rows in results
  tables_list_ordered <- reorder_results(tables_list)
  
  # Check any tables with "index" is renamed to "Gene"
  tables_list_ordered_renamed <- lapply(tables_list_ordered, function(df) {
    if ("index" %in% colnames(df)) {
      df <- rename(df, "Gene" = "index")
    }
    return(df)
  })
  
  # Create .cols sidecar meta data file
  Map(function(path, table) {
    # Create .cols sidecar meta data file
    create_cols_meta(
      file_name = path,
      table_variable_name = table,
      colname_descriptions = colname_descriptions
    )
  }, new_paths, tables_list_ordered_renamed)

  # Save the tables with the new names and p corrected values
  Map(function(df, path) {
    fwrite(df, path)
  }, tables_list_ordered_renamed, new_paths)

  # Paste to console sentence to include in manuscript
  paste_sentences(tables_list_ordered)
  
}

main(paths)
