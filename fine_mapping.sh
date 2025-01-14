# Test SuSiEx script before adding to nextflow pipeline:
# ----------------------------------------------------
# Make LD reference files (submit as jobs):
for cluster in EUR AFR SAS; do
    qsub -v cluster=$cluster -t 1-22 make-bed_hg19.sh
done

# ----------------------------------------------------
# Load interactive session for the rest of this script
qlogin -l h_vmem=32G -pe sharedmem 8

cd /exports/eddie/scratch/aedmond3/GitRepos/antidep-gwas

export PATH=$PATH:/exports/igmm/eddie/GenScotDepression/amelia/packages/plink2
module load roslin/bcftools/1.20
module load roslin/samtools/1.9

# ----------------------------------------------------
# ----- Calculate effective sample size for each meta-analysis sumstats
# Edit the meta.nf to include this as a process:
R

library(dplyr)

sumstats <- read.csv("sumstats.csv")
sumstats %>% 
  mutate(neff = 4/(1/cases + 1/controls)) %>%
  group_by(pheno, cluster) %>%
  summarise(neff = sum(neff))

write.csv(sumstats, "sumstats.Neff.csv", quote = F, row.names = F)

quit()

# ----------------------------------------------------
# ----- Process sumstats -----------------------------
# Convert sustats into correct format (convert to nextflow process)
for cluster in EUR AFR SAS; do
echo -e "SNP\tA1\tA2\tfreq\tBETA\tSE\tP\tN\tCHR\tBP\tNE" > test/fixed-N06A-${cluster}.human_g1k_v37.neff08.txt
bcftools query -f "%ID\\t%ALT\\t%REF\\t[%AFCON]\\t[%ES]\\t[%SE]\\t[%LP]\\t[%SS]\\t%CHROM\\t%POS\\t[%NE]" liftover/fixed-N06A-${cluster}.human_g1k_v37.vcf.gz | awk -v OFS='\t' -v neff_threshold=0.8 '$11 >= neff_threshold * $11 {print $1, $2, $3, $4, $5, $6, 10^-($7), $8, $9, $10, $11}' >> test/fixed-N06A-${cluster}.human_g1k_v37.neff08.txt
done

# Count the number of times zero occurs in the standard error column
# These rows where SE == 0 should be removed for SuSiEx to run
# Please see: https://github.com/getian107/SuSiEx/issues/20

# This uses awk to count where column number 6 (which is the SE column) matches "0"
# count++ increments the count for each occurrence of 0.
awk '$6 == 0 { count++ } END { print count }' test/fixed-N06A-EUR.human_g1k_v37.neff08.txt
awk '$6 == 0 { count++ } END { print count }' test/fixed-N06A-AFR.human_g1k_v37.neff08.txt
awk '$6 == 0 { count++ } END { print count }' test/fixed-N06A-SAS.human_g1k_v37.neff08.txt

# Note - check why there are occurances of a zero SE to begin with
# Send 2-3 rsIDs to Mark to trace back to find where these zeros are introduced
awk '$6 == 0 { print; if (++count == 5) exit }' test/fixed-N06A-EUR.human_g1k_v37.neff08.txt
# rs6703008, rs3128108
awk '$6 == 0 { print; if (++count == 5) exit }' test/fixed-N06A-AFR.human_g1k_v37.neff08.txt
# rs7526076, rs12132100
awk '$6 == 0 { print; if (++count == 5) exit }' test/fixed-N06A-SAS.human_g1k_v37.neff08.txt
# rs880050, rs2180311

SNPs="rs6703008 rs3128108 rs7526076 rs12132100 rs880050 rs2180311"

for SNP in $SNPs; do
zcat liftover/fixed-N06A-EUR.human_g1k_v37.vcf.gz | grep $SNP
zcat meta/fixed-N06A-EUR.meta.gz | grep $SNP
zcat sumstats/*N06A-EUR*.gz | grep $SNP
done

# Remove rows where SE == 0 in R:
# Also create a column for CPID
# ----
R

library(data.table)
library(dplyr)

lapply(c("EUR", "AFR", "SAS"), function(cluster){
  # Read in sumsstats and remove rows where SE = 0, create new col for CPID
  sumstats <- fread(paste0("test/fixed-N06A-", cluster, ".human_g1k_v37.neff08.txt")) %>% 
    filter(SE != 0) %>%
    mutate(CPID = str_glue("{CHR}:{BP}:{A2}:{A1}"))
  # Rewrite sumstats with suffix "noZero"
  fwrite(sumstats, paste0("test/fixed-N06A-", cluster, ".human_g1k_v37.neff08_noZero.txt"), sep = "\t")
  # Return number of rows in sumstats
  return(nrow(sumstats))
})

quit()
# ----

# ----------------------------------------------------
# ----- Identify regions to perform fine mapping on --
# Follow the same methods as in the SuSiEx paper:

# "in cross-population fine-mapping, we analyzed loci that reached genome-wide significance in at least one of the population-specific GWAS "

# "Loci definition.
# We used a 6-way LD clumping-based method to define genomic loci, using 1KG data as the LD reference for clumping. CEU, GBR, TSI, FIN and IBS were combined as the reference for the EUR population; ESN, GWD, LWK, MSL and YRI were combined as the reference for the AFR population; CHB, CHS, CDX, JPT and KHV were combined as the reference for the EAS population. We extracted all variants with M⁢A⁢F  >0.5%, and for each of the 25 traits, performed the LD clumping in the three populations using the corresponding reference panel and PLINK. 
# To include loci that reached genome-wide significance (𝑃  <5⁢E −8) only in the meta-analysis, we further performed clumping for the meta-GWAS across the three populations, using the three reference panels, respectively. For each clumping, we set the p-value threshold of the leading variant as 5E-8 (--clump-p1) and the threshold of the tagging variant as 0.05 (--clump-p2), and set the LD threshold as 0.1 (--clump-r2) and the distance threshold as 250 kb (--clump-kb). 
# We then took the union of the 6-way LD clumping results and extended the boundary of each merged region by 100 kb upstream and downstream. 
# Finally, we merged adjacent loci if the LD (𝑟2) between the leading variants was larger than 0.6 in any LD reference panel. "

# The changes we made to this method: using a larger clumping window of 1Mb and then not doing the merging of loci if r2 > 0.6

# --------------------
# We will perform clumping on each ancestry separately
# Let's start with the EUR ancestry as this is the only ancestry that had GWS SNPs in our GWAS, and our aim is to find the BP windows around the GWS SNPs we want to perform SuSiEx on.

# I'm also going to perform clumping on the SNPs AFTER SE == 0 rows are removed

plink2 \
  --clump test/fixed-N06A-EUR.human_g1k_v37.neff08_noZero.txt \
  --clump-id-field SNP \
  --clump-p-field P \
  --clump-p1 5e-8 \
  --clump-p2 0.05 \
  --clump-r2 0.1 \
  --clump-kb 1000 \
  --pgen reference/ukb_imp_v3.qc.geno02_EUR.pgen \
  --psam reference/ukb_imp_v3.qc.geno02_EUR.psam \
  --pvar reference/ukb_imp_v3.qc.geno02_EUR.pvar.zst \
  --out test/fixed-N06A-EUR.human_g1k_v37.neff08_noZero \
  --threads 8

# fixed-N06A-EUR.human_g1k_v37.neff08_noZero.log:
# Warning: 5001 top variant IDs in --clump file missing from main dataset.  IDs
# written to test/fixed-N06A-EUR.human_g1k_v37.neff08_noZero.clumps.missing_id .
# --clump: 1479 clumps formed from 29395 index candidates.  


# --------------------
# Load the clumping results into R to figure out the min and max BP for each clump region following a similar method as in SuSiEx
# ----
R

library(data.table)
library(dplyr)

# Read in clumps results
clumped_data <- fread("test/fixed-N06A-EUR.human_g1k_v37.neff08_noZero.clumps")

# How many clumped regions do we have?
nrow(clumped_data)
# 63

loci <- clumped_data %>% 
  dplyr::select(CHR= `#CHROM`, BP_START=POS, SNP=ID)

# Extend each locus by 100 kb upstream and downstream
loci <- loci %>%
  mutate(BP_END = BP_START + 100000,  # 100 kb downstream
         BP_START = pmax(0, BP_START - 100000)) %>% # 100 kb upstream (ensure it doesn't go negative) 
         arrange(CHR, BP_START) %>%
  relocate(CHR, BP_START, BP_END, SNP)

# Save lead SNPs as a list to calculate LD between them
write.table(loci$SNP, "test/fixed-N06A-EUR.human_g1k_v37.neff08_noZero.leading_snps",
row.names = F, quote = F, sep = "\t", col.names = F)

# Save the table of min and max BP positions for SuSiEx
write.table(loci, "test/fixed-N06A-EUR.human_g1k_v37.neff08_noZero.clumpRanges",
row.names = F, quote = F, sep = "\t")

quit()

# ----------------------------------------------------
# ----- Run SuSiEX (convert to nextflow process) -----

# BP ranges are from test/fixed-N06A-EUR.human_g1k_v37.neff08_noZero.clumpRanges
# sample sizes are taken from sumstats.Neff.csv

# Unsure whether to run these as separate jobs, as that's a lot of jobs
# Or loop them within a job
# It is quite quick to run SuSiEx so i think looping within one job will be better

# ----------------------------------------------------
# File containing CHR, BP_START, BP_END
CLUMP_RANGES_FILE="test/fixed-N06A-EUR.human_g1k_v37.neff08_noZero.clumpRanges"

start=`date +%s`

tail -n +2 "$CLUMP_RANGES_FILE" | while IFS=$'\t' read -r line; do
    # Get CHR, BP_START, BP_END by using awk to match column names indices
    # Really important we put these cols in the correct order in the R scripts making the "loci" object
    CHR=$(echo "$line" | awk '{print $1+0}')  
    BP_START=$(echo "$line" | awk '{print $2+0}')
    BP_END=$(echo "$line" | awk '{print $3+0}')
    
    echo "Processing CHR: $CHR, BP_START: $BP_START, BP_END: $BP_END"

    ../SuSiEx/bin/SuSiEx \
      --sst_file=test/fixed-N06A-EUR.human_g1k_v37.neff08_noZero.txt,test/fixed-N06A-AFR.human_g1k_v37.neff08_noZero.txt,test/fixed-N06A-SAS.human_g1k_v37.neff08_noZero.txt \
      --n_gwas=667771,39866,5814 \
      --ref_file=reference/ukb_imp_v3.qc.geno02.mind02_EUR_${CHR},reference/ukb_imp_v3.qc.geno02.mind02_AFR_${CHR},reference/ukb_imp_v3.qc.geno02.mind02_SAS_${CHR} \
      --ld_file=fineMapping/ld/EUR_${CHR},fineMapping/ld/AFR_${CHR},fineMapping/ld/SAS_${CHR} \
      --out_dir=./fineMapping/results \
      --out_name=SuSiEx.EUR.AFR.SAS.output.cs95_${CHR}:${BP_START}:${BP_END} \
      --level=0.95 \
      --pval_thresh=1e-5 \
      --chr=$CHR \
      --bp=$BP_START,$BP_END \
      --maf=0.005 \
      --snp_col=1,1,1 \
      --chr_col=9,9,9 \
      --bp_col=10,10,10 \
      --a1_col=2,2,2 \
      --a2_col=3,3,3 \
      --eff_col=5,5,5 \
      --se_col=6,6,6 \
      --pval_col=7,7,7 \
      --mult-step=True \
      --plink=../SuSiEx/utilities/plink \
      --keep-ambig=True |& tee fineMapping/logs/SuSiEx.EUR.AFR.SAS.output.cs95_${CHR}:${BP_START}:${BP_END}.log

done < "$CLUMP_RANGES_FILE"
end=`date +%s`

runtime=$((end-start))
echo $runtime


# ----------------------------------------------------
# ----- Explore results ---------------------------------

# A .summary file, a .cs file and a .snp file will be written to the specified output directory.

# If the varitional algorithm did not converge, "FAIL" will be written to both the .summary
# file and the .cs file. If no credible set was identified at the specified coverage level 
# after purity and marginal p-value filtering, "NULL" will be written to both the .summary 
# file and the .cs file.

# Check we have the same number of output files as the number of BP ranges
ls fineMapping/results/*.summary | wc -l 

# Check which files contain NULL
grep -o "NULL" fineMapping/results/*.summary

# Count how many files contain NULL
grep -o "NULL" fineMapping/results/*.summary | wc -l 
# 46 files contain the word null

# Check which files contain FAIL
grep -l "FAIL" fineMapping/results/*.summary 
# no files

# Check which files contain results
grep -L "NULL" fineMapping/results/*.summary | wc -l
# 19 files with results


# From the combined .summary files:
# Load them into R and format "#CHR:BP:BP" into new separate columns 
# Explore the following:
# - How many credible sets per fine mapped region
# - How does the max PIP vary between top SNPs
# - Plot the p-val of the top SNPs

# With the .cs and .snp files:
# - Location of SNPs in credible sets (locus zoom plot - Mark do you have code for this?)
#   - May also need data from .snp files (contains information for all the SNPs that are used in the fine-mapping algorithm)
# - Plot just the SNPs in the credible sets to see how close their PIP is, this will help determine the confidence in the causal SNP.

# ------------------------------------------
R

library(cowplot)
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)
library(stringr)
# -------------
# Define a function to process each file
process_file <- function(file, file_ext) {

  # Error handling using tryCatch to handle any potential file read errors
  tryCatch({

    # Process .summary files
    if (file_ext == "summary") {
      # Read the file as lines
      results <- fread(file, skip = 1) 
      # Process .cs or .snp files
    } else if (file_ext == "cs" | file_ext == "snp") {
      # Read in the file as a data.frame
      results <- fread(file)
    }
      # Check if the file has no rows (empty or failed file), in which case we skip it
      # Some of the SNP files are not empty, and contain only CHR and BP cols
      # These can be removed too because their corresponding .cs files contain NULL
      if (nrow(results) == 0 | (file_ext == "snp" & ncol(results) == 2)) {
        message(paste("Skipping file (it has no results):", file))
        return(NULL)
      }
    
    # Extract chromosome (CHR), start position (BP_START), and end position (BP_END)
    # These values are extracted from the file name using regular expressions
    # The regex pattern matches filenames like: "cs95_1:1000:2000.cs" or "cs95_1:1000:2000.snp"
    chr_parts <- str_match(file, "(\\d+):(\\d+):(\\d+)\\.(cs|snp|summary)")
      
    # Add the new columns to the data
    results <- results %>%
      mutate(
        CHR = as.numeric(chr_parts[2]),
        BP_START = as.integer(chr_parts[3]),
        BP_END = as.integer(chr_parts[4])
      )

    return(results)

  }, error = function(e) {
    message(paste("Error processing file:", file, "Error:", e$message))
    return(NULL)
  })
}

format_results <- function(path){
  file_extensions <- c("summary", "cs", "snp")

  results <- lapply(file_extensions, function(ext){
            paths <- list.files(path, 
              pattern = paste0("\\.", ext ,"$"), 
              full.names = TRUE)

            processed_files <- lapply(paths, function(path){
              process_file(file = path, file_ext = ext)
              })
            # Clean up NULL values (failed files)
            processed_files <- processed_files[!sapply(processed_files, is.null)]

            return(processed_files)
            })

  names(results) <- file_extensions
  results$summary <- do.call(rbind, results$summary)
  return(results)
}

path_to_susiex_results <- "fineMapping/results/"
results <- format_results(path_to_susiex_results)

# Check that each processed file type has the same number of fine mapped regions 
nrow(results$summary) == length(results$cs) && length(results$cs) == length(results$snp)


##############################
##### PLOTS ##################
# ----------------
# ------- Plot the relationship between:
# CS_LENGTH - number of SNPs in the credible set
# CS_PURITY - purity of the credible set
# MAX_PIP - Maximum posterior inclusion probability (PIP) in the credible set.
# ----------------

plotLengthPurityMaxPIP <-  function(summary_results = results$summary){
  summary_results %>%
  arrange(as.numeric(CHR)) %>%
  mutate(CHR = factor(CHR, levels = unique(CHR))) %>%
  ggplot(., aes(x = CS_PURITY, y = MAX_PIP, size = CS_LENGTH)) +
    geom_point(alpha = 0.6, aes(colour = CHR)) +  
    scale_size_continuous(range = c(1, 10)) +  
    labs(
      x = "Credible Set Purity",
      y = "Maximum Posterior Inclusion Probability",
      size = "Number of SNPs\nin credible set\n(min=1; max=170)",
      colour = "Chromosome"
    ) +
    theme_minimal() +
    theme(
      # Arrange the legends side by side and wrap CHR legend into two columns
      legend.position = "right", # Place legends at the bottom
      legend.box = "horizontal",  # Arrange legends horizontally
    ) +
    guides(
      # Set the legend for CHR (color) to have 2 columns
      colour = guide_legend(ncol = 2, order = 2), 
      size = guide_legend(order = 1),  # Ensure size legend appears first
    )
}

png("fineMapping/plots/length_purity_maxPIP.png", width = 1000, height = 600, res = 150)
  plotLengthPurityMaxPIP(results$summary)
dev.off()

# ----------------
# -------- Plot the probability the top SNP in the credible set is causal in each ancestry 
# ----------------
# Determine which ancestries were used in the fine mapping,
# do this by seeing where NA occurs in eg. -LOG10P col
# maybe also combine with POST-HOC_PROB_POP1,2,3 cols:
# POST-HOC_PROB_POP${i}: Post hoc probability credible set manifest causal in population ${i}.

ancestries <- c("EUR", "AFR", "SAS") # this variable will be definied further up in the nextflow pipeline, as input for susiex

plotAncestryCausal <-  function(summary_results = results$summary){
  summary_results %>%
    # Create a new column called ANCESTRY which tells us which ancestries have data on that SNP
    # Many of the columns containing info from all ancestries
    # are separated by commas for each ancestry, in the order of the susiex command
    mutate(ANCESTRY = ALT_ALLELE) %>%
    separate(ANCESTRY, into = ancestries, sep = ",") %>%
    # If there is no data for that ancestry there is an NA
    mutate(ANCESTRY = pmap_chr(across(all_of(ancestries)), function(...) {
      values <- c(...)
      present_ancestries <- ancestries[values != "NA"]  # Check which are not NA
      if (length(present_ancestries) > 0) {
        paste(present_ancestries, collapse = ",")  # Concatenate present ancestries
      } else {
        "None"  # If all are NA (this should not happen)
      }
    })) %>%
    # remove the indiviudal ancestry columnms now we have the column we want
    dplyr::select(-c("EUR", "AFR", "SAS")) %>%
    # create a column for location on the genome
    mutate(LOCATION = str_glue("{CHR}:{BP_START}:{BP_END}")) %>%
    # Ensure CHR and BP are in correct order on x-axis
    arrange(as.numeric(CHR), as.numeric(BP_START)) %>%
    mutate(LOCATION = factor(LOCATION, levels = LOCATION)) %>%
    # Reshape the data so that POST-HOC_PROB_POP columns are turned into long format
    pivot_longer(
      cols = starts_with("POST-HOC_PROB_POP"),  
      names_to = "POST_HOC_PROB_POP_ANCESTRY",  
      values_to = "POST_HOC_PROB_POP_CS"        
    ) %>%
    mutate(POST_HOC_PROB_POP_ANCESTRY = ancestries[as.integer(str_extract(POST_HOC_PROB_POP_ANCESTRY, "\\d+"))]) %>%
    mutate(POST_HOC_PROB_POP_ANCESTRY = factor(POST_HOC_PROB_POP_ANCESTRY, levels = unique(POST_HOC_PROB_POP_ANCESTRY))) %>%
    ggplot(data = .) +
      geom_point(aes(x = LOCATION, 
                    y = POST_HOC_PROB_POP_CS, 
                    shape = POST_HOC_PROB_POP_ANCESTRY,
                    color = ANCESTRY)) + 
      labs(
        x = "Fine mapped region",
        y = "Post-hoc probability credible set is causal",
        shape = "Ancestry",
        color = "Ancestries of LD references\nused in SuSiEx"
      ) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}


png("fineMapping/plots/POST-HOC_PROB_POP.png", width = 1000, height = 800, res = 150)
  plotAncestryCausal(results$summary)
dev.off()

####################################
### LOCUS ZOOM PLOTS ###############

# Loop over fine mapped regions, creating a new plot file for each containing:

# Plot a region plot of -log10(p) vs SNP for each ancestry (1 ancestry plot per row)
# Another option is to overlay ancestries on top of each other and not colour by R2, but that might be messy?
# Bottom row is PIP vs SNP plot (the PIP is the overall PIP from susieX, a combination of all ancestries)s

# Define fine mapped regions:
# These will be each element of results$snp and results$cs
# They should be in the same order

# ------
# To get the p-values for all the SNPs we need to load in the gwas sumstats
# Note not all these SNPs were included in the fine mapping, 
# perhaps indicate the SNPs that were included on the plot
# test/fixed-N06A-EUR.human_g1k_v37.neff08_noZero.txt,test/fixed-N06A-AFR.human_g1k_v37.neff08_noZero.txt,test/fixed-N06A-SAS.human_g1k_v37.neff08_noZero.txt
sumstats <- lapply(ancestries, function(ancestry){ # change input to actual path in nextflow process
  fread(paste0("test/fixed-N06A-",ancestry,".human_g1k_v37.neff08_noZero.txt")) 
  })
names(sumstats) <- ancestries

# Code of for loop too long so put it in a separate file
source("fine_mapping_plots.R")

# Notes:
# topr adds it's own gene annotations.. can we override this with the mBAT-combo results?