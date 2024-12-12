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

# Creates separate bim, bed, fam files for each ancestry from the UK Biobank pgen files (hg19)
# geno 0.02 = filters out variants with missing genotype data
# mind 0.02 = filters out individuals with missing genotype data
# keep-col-match keeps IDs where the string matches the ancestry cluster (ie. EUR, AFR or SAS)

cd /exports/eddie/scratch/$USER/GitRepos/antidep-gwas

# Create a file for each ancestry in UKB (EUR AFR SAS) with columns for UKB ID (FID and IID)
# Too many individuals were dropped when mind was 0.02, so try filtering geno first then mind
/gpfs/igmmfs01/eddie/GenScotDepression/amelia/packages/plink2 \
--pfile /exports/igmm/eddie/GenScotDepression/data/ukb/genetics/impv3_pgen/ukb_imp_v3.qc \
--keep-col-match reference/ukb-ld-ref_ancestry.id $cluster \
--geno 0.02 \
--make-pgen 'vzs' \
--out reference/ukb_imp_v3.qc.geno02_${cluster} \
--threads 4

/gpfs/igmmfs01/eddie/GenScotDepression/amelia/packages/plink2 \
--pfile reference/ukb_imp_v3.qc.geno02_${cluster} 'vzs' \
--mind 0.02 \
--make-bed \
--out reference/ukb_imp_v3.qc.geno02.mind02_${cluster} \
--threads 4


