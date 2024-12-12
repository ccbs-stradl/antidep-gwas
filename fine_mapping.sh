# Test SuSiEx script before adding to nextflow pipeline:
# ----------------------------------------------------
qlogin -l h_vmem=32G

cd /exports/eddie/scratch/aedmond3/GitRepos/antidep-gwas

module load roslin/bcftools/1.20
module load roslin/samtools/1.9

# ----------------------------------------------------
# ----- Process sumstats -----------------------------
# Convert sustats into correct format (convert to nextflow process)
for cluster in EUR AFR SAS; do
echo -e "SNP\tA1\tA2\tfreq\tBETA\tSE\tP\tN\tCHR\tBP\tNE" > test/${cluster}.sumstats.txt
bcftools query -f "%ID\\t%ALT\\t%REF\\t[%AFCON]\\t[%ES]\\t[%SE]\\t[%LP]\\t[%SS]\\t%CHROM\\t%POS\\t[%NE]" liftover/fixed-N06A-${cluster}.human_g1k_v37.vcf.gz | awk -v OFS='\t' -v neff_threshold=0.8 '$11 >= neff_threshold * $11 {print $1, $2, $3, $4, $5, $6, 10^-($7), $8, $9, $10, $11}' >> test/${cluster}.sumstats.txt
done

# ----------------------------------------------------
# ----- Identify lead SNPs for each ancestry ---------
# "in cross-population fine-mapping, we analyzed loci that reached genome-wide significance in at least one of the population-specific GWAS "

# First count how many genome-wide significant (GWS) SNPs are in each ancestry
for cluster in EUR AFR SAS; do
cat test/${cluster}* | awk '$7 <= 5e-8' | wc -l
done
# EUR: 2151; AFR: 0; SAS:0

# Create Supplementary Table of all 2,151 GWS SNPs with the following columns:
# SNP, EUR_P_VAL, AFR_P_VAL, SAS_P_VAL
cat test/* \
| awk '$7 <= 5e-8' \
| awk 'a[$1]++ == 0' \
| awk '{print $1}' \
| sort > fineMapping/GWS_SNPs_tmp.txt

for cluster in EUR AFR SAS; do
cat test/${cluster}* \
| awk -v cluster="${cluster}" '$1 != "." {print $1, $7}' \
| sort > fineMapping/GWS_SNPs_tmp_${cluster}.txt
done

echo -e "SNP\tEUR_P_VAL\tAFR_P_VAL\tSAS_P_VAL" > fineMapping/GWS_SNPs_PVals.txt
join -a1 -e "NA" fineMapping/GWS_SNPs_tmp.txt fineMapping/GWS_SNPs_tmp_EUR.txt \
| join -a1 -e "NA" - fineMapping/GWS_SNPs_tmp_AFR.txt \
| join -a1 -e "NA" - fineMapping/GWS_SNPs_tmp_SAS.txt >> fineMapping/GWS_SNPs_PVals.txt

# join: fineMapping/GWS_SNPs_tmp_EUR.txt:16: is not sorted: rs1000002 0.9628  
# join: fineMapping/GWS_SNPs_tmp_AFR.txt:21: is not sorted: rs1000002 0.785001
# join: fineMapping/GWS_SNPs_tmp_SAS.txt:17: is not sorted: rs1000002 0.7926  
# join: input is not in sorted order
# join: input is not in sorted order
# join: input is not in sorted order

# Not convinced this has worked correctly, tempted to read it into R or Python and wrangle there
head fineMapping/GWS_SNPs_PVals.txt


# We also need to save the chromosome and base position for each GWS SNP
# then calculate the min and max base position value of the region for fine mapping
# ie. define the region we want to fine map around lead SNPs
# "1-Mb regions" used in SuSiEx paper simulations, ie. BP ± 500,000 for min/max values
# what if BP ±500,000 is outside the range of the chromosome?
# Also there's a lot of SNPs that overlap...
# Unsure how best to get the BP start and end values for each fine mapping region
# Let's do this in R.
R

library(data.table)
library(dplyr)
library(purrr)

# Read in GWS SNPs
GWS_SNPs <- fread("fineMapping/GWS_SNPs_tmp.txt", header = FALSE) %>% pull(1)

# Read in sumstats and subset to GWS SNPs
sumstats <- lapply(c("EUR", "AFR", "SAS"), function(cluster){
  fread(paste0("test/", cluster, ".sumstats.txt")) %>%
  filter(SNP %in% GWS_SNPs) %>%
  rename_with(~ paste0(cluster, "_", .), -SNP)  # Add cluster prefix to column names, except for SNP
})

# Check if they all contain those SNPs
lapply(sumstats, nrow)
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

quit()

rm fineMapping/*_tmp_*

# ----------------------------------------------------
# ----- Calculate LD matrices for each region that will be fine mapped
# Use plink2 and use the SuSiEx script as a guide: https://github.com/getian107/SuSiEx/blob/master/utilities/SuSiEx_LD.py

# Extract list of SNPs for region of interest, using CHR and BP parameters in the .bim file of reference/ukb_imp_v3.qc_EUR. Refered to below as $CHR and reference/ukb_imp_v3.qc_EUR.snp
# or we could use  and --to-bp below instead

# Create bed files for that region
CHR=1
BP_START=7314654
BP_END=8314677
cluster=EUR

../SuSiEx/utilities/plink \
--bfile reference/ukb_imp_v3.qc_${cluster} \
--keep-allele-order \
--chr $CHR \
--from-bp $BP_START \
--to-bp $BP_END \
--make-bed \
--out reference/ukb_imp_v3.qc_${cluster}_ref

# Create LD matrix (if changing to PLINK2 change to "--r-unphased" instead of "--r")
../SuSiEx/utilities/plink \
--bfile reference/ukb_imp_v3.qc_${cluster}_ref \
 --keep-allele-order \
 --r square bin4 \
 --out reference/ukb_imp_v3.qc_${cluster}

# Calculate allele frequencies:
../SuSiEx/utilities/plink \
--bfile reference/ukb_imp_v3.qc_${cluster}_ref \
--keep-allele-order \
--freq \
--out reference/ukb_imp_v3.qc_${cluster}


# Loop over all ancestries:
#reference/ukb_imp_v3.qc_EUR, reference/ukb_imp_v3.qc_AFR,reference/ukb_imp_v3.qc_SAS


# Calculate and save how many variants are in each region (for Supplementary Table)





# ----------------------------------------------------
# ----- Run SuSiEX (convert to nextflow process) -----
../SuSiEx/bin/SuSiEx \
  --sst_file=test/EUR.sumstats.txt,test/AFR.sumstats.txt,test/SAS.sumstats.txt \
  --n_gwas=5800000,48000,7300 \
  --ld_file=reference/ukb_imp_v3.qc_EUR_LD,reference/ukb_imp_v3.qc_AFR_LD,reference/ukb_imp_v3.qc_SAS_LD \
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
  --keep-ambig=True |& tee fineMapping/SuSiEx.EUR.AFR.SAS.output.cs95.log

../SuSiEx/bin/SuSiEx \
  --sst_file=test/EUR.sumstats.txt,test/SAS.sumstats.txt \
  --n_gwas=5800000,7300 \
  --ref_file=reference/ukb_imp_v3.qc_EUR,reference/ukb_imp_v3.qc_SAS \
  --ld_file=fineMapping/EUR,fineMapping/SAS \
  --out_dir=./fineMapping \
  --out_name=SuSiEx.EUR.SAS.output.cs95 \
  --level=0.95 \
  --pval_thresh=1e-5 \
  --chr=$CHR \
  --bp=$BP_START,$BP_END \
  --maf=0.005 \
  --snp_col=1,1 \
  --chr_col=9,9 \
  --bp_col=10,10 \
  --a1_col=2,2 \
  --a2_col=3,3 \
  --eff_col=5,5 \
  --se_col=6,6 \
  --pval_col=7,7 \
  --plink=../SuSiEx/utilities/plink \
  --mult-step=True \
  --keep-ambig=True |& tee fineMapping/SuSiEx.EUR.SAS.output.cs95.log


# Error when using plink2 rather than plink
../SuSiEx/bin/SuSiEx \
  --sst_file=test/EUR.sumstats.txt,test/AFR.sumstats.txt,test/SAS.sumstats.txt \
  --n_gwas=5800000,48000,7300 \
  --ref_file=reference/ukb_imp_v3.qc_EUR,reference/ukb_imp_v3.qc_AFR,reference/ukb_imp_v3.qc_SAS \
  --ld_file=fineMapping/EUR,fineMapping/AFR,fineMapping/SAS \
  --out_dir=./fineMapping \
  --out_name=SuSiEx.EUR.AFR.SAS.output.cs95 \
  --level=0.95 \
  --pval_thresh=1e-5 \
  --chr=1 \
  --bp=7314654,8314677 \
  --maf=0.005 \
  --snp_col=1,1,1 \
  --chr_col=9,9,9 \
  --bp_col=10,10,10 \
  --a1_col=2,2,2 \
  --a2_col=3,3,3 \
  --eff_col=5,5,5 \
  --se_col=6,6,6 \
  --pval_col=7,7,7 \
  --plink=/gpfs/igmmfs01/eddie/GenScotDepression/amelia/packages/plink2 \
  --mult-step=True \
  --keep-ambig=True |& tee fineMapping/SuSiEx.EUR.AFR.SAS.output.cs95_PLINK_DEBUG.log

 # The version of plink used in SuSiEx requires "--r" to be available in plink. I think this has been removed as an option in plink2. As using plink2 (as above) gives the error:
# "Error: Unrecognized flag ('--r')."


# ----------------------------------------------------
# With the above code there are a few things to change:
# 1. This is example code so "chr" and "bp" need changing to the regions needing fine mapping.
# 2. But more importantly, the LD matrices do not seem to be being calculated for SAS. (See log file). 
  # Solution to this is to use PLINK to calculate LD matrices directly (using the genome wide signifcant SNP regions) and parse them to "--ld-file", if doing this remove flags for "--ref-file" and "--plink".
  # See their script to calculate LD matrices here: https://github.com/getian107/SuSiEx/blob/master/utilities/SuSiEx_LD.py, maybe adapt it to use plink2.
# 3. Also need to check what to use as sample size for "--n_gwas", is this the sample size stated in the abstract of these studies or is this the effective sample size for the study?

