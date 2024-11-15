#$ -N make-bed_hg19
#$ -l h_vmem=32G
#$ -l rl9=true
#$ -pe sharedmem 4
#$ -l h_rt=8:00:00
#$ -e logs/
#$ -o logs/
#$ -cwd
#$ -m beas
#$ -M aedmond3@ed.ac.uk

# Use UKBiobank genetic data for creating LD references as these are in the hg19 format
# Subset pgen files to specific ancestries

cd /exports/eddie/scratch/$USER/GitRepos/antidep-gwas

# Create a file for each ancestry in UKB (EUR AFR SAS) with columns for UKB ID (FID and IID)
for cluster in EUR AFR SAS; do
/gpfs/igmmfs01/eddie/GenScotDepression/amelia/packages/plink2 \
--pfile /exports/igmm/eddie/GenScotDepression/data/ukb/genetics/impv3_pgen/ukb_imp_v3.qc \
--keep /exports/igmm/eddie/GenScotDepression/data/ukb/genetics/input_filters/similarity_clusters/ld_ref/ukb-ld-ref_${cluster}.id \
--maf 0.01 \
--mac 100 \
--geno 0.02 \
--mind 0.02 \
--make-bed \
--out reference/ukb_imp_v3.qc_${cluster} \
--threads 4
done


