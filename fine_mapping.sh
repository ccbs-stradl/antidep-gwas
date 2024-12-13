# Test SuSiEx script before adding to nextflow pipeline:
# ----------------------------------------------------
# Make reference files:
for cluster in EUR AFR SAS; do
  qsub -v cluster=$cluster make-bed_hg19.sh
done

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
qlogin -l h_vmem=32G

cd /exports/eddie/scratch/aedmond3/GitRepos/antidep-gwas

module load roslin/bcftools/1.20
module load roslin/samtools/1.9

# ----------------------------------------------------
# ----- Process sumstats -----------------------------
# Convert sustats into correct format (convert to nextflow process)
for cluster in EUR AFR SAS; do
echo -e "SNP\tA1\tA2\tfreq\tBETA\tSE\tP\tN\tCHR\tBP\tNE" > test/fixed-N06A-${cluster}.human_g1k_v37.neff08.txt
bcftools query -f "%ID\\t%ALT\\t%REF\\t[%AFCON]\\t[%ES]\\t[%SE]\\t[%LP]\\t[%SS]\\t%CHROM\\t%POS\\t[%NE]" liftover/fixed-N06A-${cluster}.human_g1k_v37.vcf.gz | awk -v OFS='\t' -v neff_threshold=0.8 '$11 >= neff_threshold * $11 {print $1, $2, $3, $4, $5, $6, 10^-($7), $8, $9, $10, $11}' >> test/fixed-N06A-${cluster}.human_g1k_v37.neff08.txt
done

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

# ---------------
# Check if GWS SNP fall within .clumped.ranges output from meta.nf
# because our gwas sumstats are now on hg19 whereas the output in meta is in Gr38.
# I'm not sure if this would mean the clumped ranges are different or not

clumped_ranges <- fread("meta/fixed-N06A-EUR.clumped.ranges")
sort(clumped_ranges$POS)
# gives 110 ranges, we need to check that all our GWS SNPs fall within all these ranges
# even by eye i can see that there are GWS SNPs on CHR 1 but no CHR 1 in clumped ranges



quit()
# ----

rm fineMapping/*_tmp_*


# ----------------------------------------------------
# ----- Run SuSiEX (convert to nextflow process) -----

CHR=1
BP_START=7314654
BP_END=8314677
# use meta/*.clumped.ranges to get start and end positions and chr
# make sure they are on the correct hg19 build though
# only clumps for EUR, nothing for SAS or AFR because there weren't any GWS SNPs
# sample sizes are taken from sumstats.Neff.csv

../SuSiEx/bin/SuSiEx \
  --sst_file=test/fixed-N06A-EUR.human_g1k_v37.neff08.txt,test/fixed-N06A-AFR.human_g1k_v37.neff08.txt,test/fixed-N06A-SAS.human_g1k_v37.neff08.txt \
  --n_gwas=667771,39866,5814 \
  --ref_file=reference/ukb_imp_v3.qc_EUR,reference/ukb_imp_v3.qc_AFR,reference/ukb_imp_v3.qc_SAS \
  --ld_file=fineMapping/EUR,fineMapping/AFR,fineMapping/SAS \
  --out_dir=./fineMapping \
  --out_name=SuSiEx.EUR.AFR.SAS.output.cs95 \
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
  --keep-ambig=True |& tee fineMapping/SuSiEx.EUR.AFR.SAS.output.cs95.log


# "Error: Effect size of Line 28807 is not a number (NAN)"
# in what looks like the sumstats file
sed -n '28808p' test/EUR*
# It contains a beta value with a zero
awk '$5 == 0 { count++ } END { print count }' test/EUR*
# There's a lot of zero beta values in each sum stats
# For now let's remove them to see if the rest of SuSiEx works
awk '$5 != 0' test/fixed-N06A-EUR.human_g1k_v37.neff08.txt > test/fixed-N06A-EUR.human_g1k_v37.neff08.noZero.txt
awk '$5 != 0' test/fixed-N06A-AFR.human_g1k_v37.neff08.txt > test/fixed-N06A-AFR.human_g1k_v37.neff08.noZero.txt
awk '$5 != 0' test/fixed-N06A-SAS.human_g1k_v37.neff08.txt > test/fixed-N06A-SAS.human_g1k_v37.neff08.noZero.txt



# ----------------------------------------------------
# With the above code there are a few things to change:
# 1. This is example code so "chr" and "bp" need changing to the regions needing fine mapping. Use meta/*EUR.clumped.ranges to get base positions
# 2. SuSiEx doesn't seem to like it when effect sizes are zero in the sum stats file


