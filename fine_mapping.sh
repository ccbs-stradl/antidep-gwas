# Test SuSiEx script before adding to nextflow pipeline:

qlogin -l h_vmem=32G

cd /exports/eddie/scratch/aedmond3/GitRepos/antidep-gwas

module load roslin/bcftools/1.20
module load roslin/samtools/1.9

# Convert sustats into correct format (convert to nextflow process)
for cluster in EUR AFR SAS; do
echo -e "SNP\tA1\tA2\tfreq\tBETA\tSE\tP\tN\tCHR\tBP\tNE" > test/${cluster}.sumstats.txt
bcftools query -f "%ID\\t%ALT\\t%REF\\t[%AFCON]\\t[%ES]\\t[%SE]\\t[%LP]\\t[%SS]\\t%CHROM\\t%POS\\t[%NE]" liftover/fixed-N06A-${cluster}.human_g1k_v37.vcf.gz >> test/${cluster}.sumstats.txt
done

# Run SuSiEX (convert to nextflow process)
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
  --plink=../SuSiEx/utilities/plink \
  --mult-step=True \
  --keep-ambig=True |& tee fineMapping/SuSiEx.EUR.AFR.SAS.output.cs95.log

# With the above code there are a few things to change:
# 1. This is example code so "chr" and "bp" need changing to the regions needing fine mapping.
# 2. But more importantly, the LD matrices do not seem to be being calculated for SAS. (See log file). 
  # Solution to this is to use PLINK to calculate LD matrices directly (using the genome wide signifcant SNP regions) and parse them to "--ld-file", if doing this remove flags for "--ref-file" and "--plink".
  # Unsure how to choose the regions that need fine mapping from the sumstats of different ancestries. eg. a lead SNP for EUR sumstats may not be a lead SNP for AFR. Check literature for more info on this.