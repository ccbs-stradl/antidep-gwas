#$ -N pgen_hg19
#$ -l h_vmem=32G
#$ -l rl9=true
#$ -l h_rt=8:00:00
#$ -e logs
#$ -o logs
#$ -cwd

# Use UKBiobank genetic data for creating LD references as these are in the hg19 format
# Subset pgen files to specific ancestries

cd /exports/eddie/scratch/$USER/antidep-gwas
export PATH=$PATH:/gpfs/igmmfs01/eddie/GenScotDepression/amelia/packages/plink2

# Create a file with columns for UKB ID (FID and IID) and a column for ancestry
rm reference/ukb-ld-ref_ancestry.id
for cluster in EUR AFR SAS; do
awk -v cluster="$cluster" '{print $1, $2, cluster}' /exports/igmm/eddie/GenScotDepression/data/ukb/genetics/input_filters/similarity_clusters/ld_ref/ukb-ld-ref_${cluster}.id >> reference/ukb-ld-ref_ancestry.id
done
head reference/ukb-ld-ref_ancestry.id
tail reference/ukb-ld-ref_ancestry.id

# Save first two columns (ID columns) from ukb-ld-ref_ancestry.id as keep.txt, these are the IDs we want to keep in the pfiles
cut -d ' ' -f 1-2 reference/ukb-ld-ref_ancestry.id > reference/ukb-ld-ref_keep.txt


# Create pfile containing the UKB IDs in reference/ukb-ld-ref_keep.txt
plink2 \
    --pfile /exports/igmm/eddie/GenScotDepression/data/ukb/genetics/impv3_pgen/ukb_imp_v3.qc \
    --keep reference/ukb-ld-ref_keep.txt \
    --make-pgen 'vzs' \
    --out reference/ukb_imp_v3.qc_ancestry

# Add ancestry infomation to .psam
awk 'NR==1 {next} 
     BEGIN {print "#IID\tPAT\tMAT\tSEX"} 
     {print $2 "\t0\t0\t" $3}'  reference/ukb_imp_v3.qc_ancestry.psam > tmp_sex.psam

awk 'FNR==NR {superpop[$1] = $3; pop[$1] = $3; next} 
     FNR==1 {print $0 "\tSuperPop\tPopulation"; next} 
     {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" superpop[$1] "\t" pop[$1]}' reference/ukb-ld-ref_ancestry.id tmp_sex.psam > tmp.psam

mv tmp.psam reference/ukb_imp_v3.qc_ancestry.psam    
rm tmp_sex.psam

head reference/ukb_imp_v3.qc_ancestry.psam    
# The .psam file is now in the same format as 1000 genomes Gr38 build
