# Test SuSiEx script before adding to nextflow pipeline:
# ----------------------------------------------------
# Make reference files:
for cluster in EUR AFR SAS; do
  qsub -v cluster=$cluster make-bed_hg19.sh
done

# ----------------------------------------------------
qlogin -l h_vmem=32G

cd /exports/eddie/scratch/aedmond3/GitRepos/antidep-gwas

module load roslin/bcftools/1.20
module load roslin/samtools/1.9

# ----------------------------------------------------
# Calculate effective sample size for each meta-analysis sumstats
# Edit the meta.nf to include this as a process:
R
library(dplyr)
sumstats <- read.csv("sumstats.csv")
sumstats %>% 
  mutate(neff = 4/(1/cases + 1/controls)) %>%
  group_by(pheno, cluster) %>%
  summarise(neff = sum(neff))

mutate(neff = 1/ (1/cases + 1/controls)) |> group_by(pheno, cluster) |> summarize(neff = sum(neff))

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
awk '$6 == 0 { count++ } END { print count }' test/EUR*

# Note - check why there are occurances of a zero SE to begin with

# Remove rows where SE == 0 in R:
# Also create a column for CPID as we need this later to define BP ranges
# ----
R

library(data.table)
library(dplyr)

lapply(c("EUR", "AFR", "SAS"), function(cluster){
  # Read in sumsstats and remove rows where SE = 0
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
# ----- Identify lead SNPs for each ancestry ---------
# "in cross-population fine-mapping, we analyzed loci that reached genome-wide significance in at least one of the population-specific GWAS "

# define the region we want to fine map around lead SNPs
# "1-Mb regions" used in SuSiEx paper simulations, 
# use meta/*.clumped.ranges output from meta.nf

# First count how many genome-wide significant (GWS) SNPs are in each ancestry
for cluster in EUR AFR SAS; do
cat test/fixed-N06A-${cluster}.human_g1k_v37.neff08.txt | awk '$7 <= 5e-8' | wc -l
done
# EUR: 2151; AFR: 0; SAS:0


# ----
R

library(data.table)
library(dplyr)
library(purrr)


# ---------------
# Create table of GWS SNPs that will be fine mapped with the P-values for each ancestry:

# Read in sumstats and subset to GWS SNPs
sumstats <- lapply(c("EUR", "AFR", "SAS"), function(cluster){
  fread(paste0("test/fixed-N06A-", cluster, ".human_g1k_v37.neff08.txt")) 
})


GWS_SNPs <- lapply(sumstats, function(sumstats){
  sumstats %>%
  filter(P <= 5e-8) %>%
  pull(SNP)
}) %>% unlist()


clusters <- c("EUR", "AFR", "SAS")
sumstats <- lapply(1:3, function(i){
  sumstats[[i]] %>%
  filter(SNP %in% GWS_SNPs) %>%
  rename_with(~ paste0(clusters[i], "_", .), -SNP)  # Add cluster prefix to column names, except for SNP
})

# Check if they all contain those SNPs
lapply(sumstats, nrow)
  # [[1]]
  # [1] 2151

  # [[2]]
  # [1] 1973

  # [[3]]
  # [1] 1992
# Some SNPs are missing in AFR and SAS

# Join all the sumstats on SNP col into one table with the cluster prefix before each duplicated col name
# Left join onto EUR, to keep all GWS SNPs
sumstats <- reduce(sumstats, left_join, by = "SNP")

# Check that CHR and BP are the same across all ancestries for the same SNP (QC check of build)
# ignore missing SNPs
sum(!sumstats$EUR_CHR == sumstats$AFR_CHR, na.rm = T) # should all = 0
sum(!sumstats$EUR_CHR == sumstats$SAS_CHR, na.rm = T)
sum(!sumstats$EUR_BP == sumstats$AFR_BP, na.rm = T) # should all = 0
sum(!sumstats$EUR_BP == sumstats$SAS_BP, na.rm = T)

# Subset to keep SNP, CHR and BP, EUR_P, AFR_P, SAS_P
sumstats <- sumstats %>%
  mutate(BP_START = EUR_BP - 500000) %>%
  mutate(BP_END = EUR_BP + 500000) %>%
  select(SNP, CHR = EUR_CHR, BP = EUR_BP, BP_START, BP_END, EUR_P, AFR_P, SAS_P) 

write.table(sumstats, "fineMapping/GWS_SNPs.txt", row.names = FALSE, sep = "\t", quote = FALSE)

# Noticed a few GWS_SNPs are on chromosome X, let's just count these
sumstats %>%
  filter(CHR == "X") %>%
  nrow()
  # 110

quit()
# ----

# ----------------------------------------------------
# ----- Identify regions to perform fine mapping on --
# From the methods section of the SuSiEx paper:
# "Loci definition.
# We used a 6-way LD clumping-based method to define genomic loci, using 1KG data as the LD reference for clumping. CEU, GBR, TSI, FIN and IBS were combined as the reference for the EUR population; ESN, GWD, LWK, MSL and YRI were combined as the reference for the AFR population; CHB, CHS, CDX, JPT and KHV were combined as the reference for the EAS population. We extracted all variants with M⁢A⁢F  >0.5%, and for each of the 25 traits, performed the LD clumping in the three populations using the corresponding reference panel and PLINK45. To include loci that reached genome-wide significance (𝑃  <5⁢E −8) only in the meta-analysis, we further performed clumping for the meta-GWAS across the three populations, using the three reference panels, respectively. For each clumping, we set the p-value threshold of the leading variant as 5E-8 (--clump-p1) and the threshold of the tagging variant as 0.05 (--clump-p2), and set the LD threshold as 0.1 (--clump-r2) and the distance threshold as 250 kb (--clump-kb). We then took the union of the 6-way LD clumping results and extended the boundary of each merged region by 100 kb upstream and downstream. Finally, we merged adjacent loci if the LD (𝑟2) between the leading variants was larger than 0.6 in any LD reference panel. "


# Clump the GWAS on hg19 to see if it's different to the results in meta.nf

# Request memory and ensure plink is on the execuable path
qlogin -l h_vmem=32G -pe sharedmem 8

# Ensure plink2 is in executable path (plink 2 version: 20241222)
export PATH=$PATH:/exports/igmm/eddie/GenScotDepression/amelia/packages/plink2
cd /exports/eddie/scratch/aedmond3/GitRepos/antidep-gwas

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
  --clump-kb 250 \
  --pgen reference/ukb_imp_v3.qc.geno02_EUR.pgen \
  --psam reference/ukb_imp_v3.qc.geno02_EUR.psam \
  --pvar reference/ukb_imp_v3.qc.geno02_EUR.pvar.zst \
  --out test/fixed-N06A-EUR.human_g1k_v37.neff08_noZero \
	--threads 8

# fixed-N06A-EUR.human_g1k_v37.neff08_noZero.log:
# Warning: 5001 top variant IDs in --clump file missing from main dataset.  IDs
# written to test/fixed-N06A-EUR.human_g1k_v37.neff08_noZero.clumps.missing_id .
# --clump: 1479 clumps formed from 29395 index candidates.  

# Load the clumping results into R to figure out the min and max BP for each clump region
# ----
R

library(data.table)
library(dplyr)
library(stringr)
library(IRanges) # to detect overlap in ranges

# Read in clumps results
clumps <- fread("test/fixed-N06A-EUR.human_g1k_v37.neff08_noZero.clumps")
sumstats <- fread("test/fixed-N06A-EUR.human_g1k_v37.neff08_noZero.txt")

# For each row ie. clump set take the list of SNPs in that clumped region ie. column "SP2"
# and assign them a CPID from the sumstats file
# Loop over all rows and extract the BP_START, BP_END, CHR and lead SNP position
# check lead SNP position is within BP_START and BP_END
# check there is no overlap in regions of BP_START and BP_END

clump_ranges_df <- lapply(1:nrow(clumps), function(i){
  leadSNP <- clumps %>%
              dplyr::slice(i)

  ranges <- sumstats %>%
      filter(SNP %in% unlist(str_split(clumps$SP2[i], ","))) %>%
      group_by(CHR) %>%
      summarise(
      BP_START = min(BP),
      BP_END = max(BP)
      ) %>%
      mutate(lead_SNP = leadSNP$ID,
             lead_POS = leadSNP$POS) %>%
      mutate(inRange = lead_POS > BP_START & lead_POS < BP_END)

  return(ranges)

}) %>% do.call(rbind,.)

i <- 1475


write.csv(clump_ranges_df, "test/fixed-N06A-EUR.human_g1k_v37.neff08_noZero.clumps.pos.csv",
row.names = F, quote = F)
# clump_ranges_df <- fread("test/fixed-N06A-EUR.human_g1k_v37.neff08_noZero.clumps.pos.csv")

# Create IRanges objects for the start and end positions for each CHR
lapply(unique(clump_ranges_df$CHR), function(chr){
  clump_ranges_df_CHR <- clump_ranges_df %>%
                          filter(CHR == chr)
  ranges <- IRanges(start = clump_ranges_df_CHR$BP_START, end = clump_ranges_df_CHR$BP_END)

  return(list(chr, countOverlaps(ranges, ranges)))
})


# If it returns a vector of all 1 then there are no overlaps
# If it returns a vector where any value != 1 then we need to go back and check the clumping process
# yes it does return overlapping regions....does that matter?


# Check the lead SNPs are in fineMapping/GWS_SNPs.txt (previously identified GWS SNPs)
# They should all be in there, note that we are clumping with a less stringent P value of 0.0001, check if we need to change this to GWS

# ----------------------------------------------------
# ----- Run SuSiEX (convert to nextflow process) -----

# Tidy up directory
rm fineMapping/*_tmp_*

# BP ranges are from test/fixed-N06A-EUR.human_g1k_v37.neff08_noZero.clumps.pos.csv
# sample sizes are taken from sumstats.Neff.csv


../SuSiEx/bin/SuSiEx \
  --sst_file=test/fixed-N06A-EUR.human_g1k_v37.neff08.txt,test/fixed-N06A-AFR.human_g1k_v37.neff08.txt,test/fixed-N06A-SAS.human_g1k_v37.neff08.txt \
  --n_gwas=667771,39866,5814 \
  --ref_file=reference/ukb_imp_v3.qc.geno02.mind02_EUR,reference/ukb_imp_v3.qc.geno02.mind02_AFR,reference/ukb_imp_v3.qc.geno02.mind02_SAS \
  --ld_file=fineMapping/EUR,fineMapping/AFR,fineMapping/SAS \
  --out_dir=./fineMapping \
  --out_name=SuSiEx.EUR.AFR.SAS.output.cs95_${CHR}:${BP_START}:${$BP_END} \
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
  --keep-ambig=True |& tee fineMapping/logs/SuSiEx.EUR.AFR.SAS.output.cs95_${CHR}:${BP_START}:${$BP_END}.log

# ----------------------------------------------------
# With the above code there are a few things to change:
# 1. This is example code so "chr" and "bp" need changing to the regions needing fine mapping. Use meta/*EUR.clumped.ranges to get base positions


