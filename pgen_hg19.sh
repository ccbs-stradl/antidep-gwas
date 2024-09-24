#$ -N pgen_hg19
#$ -l h_vmem=64G
#$ -l rl9=true
#$ -l h_rt=8:00:00
#$ -e logs
#$ -o logs
#$ -cwd

# Use UKBiobank genetic data for creating LD references as these are in the hg19 format
# Subset pgen files to specific ancestries

cd /exports/eddie/scratch/$USER

plink2 \
	--pfile /exports/igmm/eddie/GenScotDepression/data/ukb/genetics/impv3_pgen/ukb_imp_v3.qc \
    --keep /exports/igmm/eddie/GenScotDepression/data/ukb/genetics/input_filters/similarity_clusters/ld_ref/ukb-ld-ref_${CLUSTER}.id \
    --make-pgen \
    --out antidep-gwas/reference/ukb_imp_v3.qc_${CLUSTER}
