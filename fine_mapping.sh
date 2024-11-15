# Test SuSiEx script before adding to nextflow pipeline:

qlogin -l h_vmem=32G

cd /exports/eddie/scratch/aedmond3/GitRepos/antidep-gwas

module load roslin/bcftools/1.20
module load roslin/samtools/1.9

# Convert sustats into correct format (convert to nextflow process)
for cluster in EUR AFR SAS; do
echo -e "SNP\tA1\tA2\tfreq\tBETA\tSE\tP\tN\tCHR\tBP\tNE" > test/${cluster}.sumstats.txt
bcftools query \
  -f "%ID\\t%ALT\\t%REF\\t[%AFCON]\\t[%ES]\\t[%SE]\\t[%LP]\\t[%SS]\\t%CHROM\\t%POS\\t[%NE]" liftover/fixed-N06A-${cluster}.human_g1k_v37.vcf.gz >> test/${cluster}.sumstats.txt
done

# Run SuSiEX (convert to nextflow process)
../SuSiEx/bin/SuSiEx \
  --sst_file=EUR.sumstats.txt,AFR.sumstats.txt,SAS.sumstats.txt \
  --n_gwas=50000,50000 \
  --ref_file=EUR,AFR,SAS \
  --ld_file=EUR,AFR,SAS \
  --out_dir=./ \
  --out_name=SuSiEx.EUR.AFR.SAS.output.cs95 \
  --level=0.95 \
  --pval_thresh=1e-5 \
  --maf=0.005
  --snp_col=1,1,1 \
  --chr_col=9,9,9 \
  --bp_col=10,10,10 \
  --a1_col=2,2,2 \
  --a2_col=3,3,3 \
  --eff_col=5,5,5 \
  --se_col=6,6,6 \
  --pval_col=7,7,7 \
  --plink=./gpfs/igmmfs01/eddie/GenScotDepression/amelia/packages/plink2 \
  --mult-step=True \
  --keep-ambig=True \
  --threads=16