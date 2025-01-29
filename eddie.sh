## Commands for running the workflows on Eddie

# load required modules and set up environment
module load roslin/nextflow/23.10.1
module load singularity/4.1.3
source /exports/applications/support/set_qlogin_environment.sh
export NXF_SINGULARITY_LIBRARYDIR=/exports/igmm/eddie/BioinformaticsResources/nfcore/singularity-images

###
### Format GWAS
###

# format antidep exposure sumstats
nextflow run format_gwas.nf -resume \
-work-dir /exports/eddie/scratch/${USER}/ad/work/format \
-c eddie.config

# convert hg38 sumstats to VCF
nextflow run vcf.nf -resume \
--sumstats "format/gwas/GRCh38/*.{txt,json,csv}" \
--chr 'chr#' \
--publish "vcf/gwas/GRCh38" \
--assembly "reference/Homo_sapiens_assembly38.{fasta,fasta.fai,dict}" \
--dbsnp "reference/dbsnp.v153.hg38.vcf.{gz,gz.tbi}" \
-work-dir /exports/eddie/scratch/${USER}/ad/work/vcf \
-c eddie.config

# convert hg19 sumstats to VCF
nextflow run vcf.nf -resume \
--sumstats "format/gwas/GRCh37/*.{txt,json,csv}" \
--chr '#' \
--publish "vcf/gwas/GRCh37" \
--assembly "reference/human_g1k_v37.{fasta,fasta.fai,dict}" \
--dbsnp "reference/dbsnp.v153.b37.vcf.{gz,gz.tbi}" \
-work-dir /exports/eddie/scratch/${USER}/ad/work/vcf \
-c eddie.config

# liftover hg19 VCFs to hg38
nextflow run liftover.nf -resume \
--sumstats "vcf/gwas/GRCh37/*.{vcf.gz,vcf.gz.tbi}" \
--source "reference/human_g1k_v37.{fasta,fasta.fai}" \
--destination "reference/Homo_sapiens_assembly38.{fasta,fasta.fai}" \
--chain reference/hg19ToHg38.over.chain.gz \
--publish liftover/gwas \
-work-dir /exports/eddie/scratch/${USER}/ad/work/liftover \
-c eddie.config 

# move lifted over file
mv liftover/gwas/*.Homo_sapiens_assembly38.vcf.* vcf/gwas/GRCh38
rename Homo_sapiens_assembly38.vcf vcf vcf/gwas/GRCh38/*.Homo_sapiens_assembly38.vcf.*
cp vcf/gwas/GRCh37/*.csv vcf/gwas/GRCh38/

###
### Run meta analysis
###

# run meta-analysis
nextflow run meta.nf -resume \
-work-dir /exports/eddie/scratch/${USER}/ad/work/meta \
-c eddie.config

# convert meta-analysis to VCF
nextflow run format_meta.nf -resume \
-work-dir /exports/eddie/scratch/${USER}/ad/work/meta \
-c eddie.config

nextflow run vcf.nf -resume \
--sumstats "format/meta/GRCh38/*.{txt,json,csv}" \
--chr 'chr#' \
--publish "vcf/meta/GRCh38" \
--assembly "reference/Homo_sapiens_assembly38.{fasta,fasta.fai,dict}" \
--dbsnp "reference/dbsnp.v153.hg38.vcf.{gz,gz.tbi}" \
-work-dir /exports/eddie/scratch/${USER}/ad/work/vcf \
-c eddie.config -with-trace

###
### Run downstream analyses
###

# Run mBAT-combo on build hg19/GRCh37
nextflow run genes.nf -resume \
-work-dir /exports/eddie/scratch/${USER}/ad/work_hg19 \
-c eddie.config \
--build 'hg19'

# Run mBAT-combo on build hg38/GRCh38
nextflow run genes.nf -resume \
-work-dir /exports/eddie/scratch/${USER}/ad/work_hg38 \
-c eddie.config \
--build 'hg38'

nextflow run metavcf.nf -resume \
-work-dir /exports/eddie/scratch/${USER}/ad/work \
-c eddie.config \
--sumstats "meta/fixed-N06A-*.meta.gz"

nextflow run liftover.nf -resume \
-work-dir /exports/eddie/scratch/${USER}/ad/work \
-c eddie.config \
--sumstats "metavcf/*.{vcf.gz,vcf.gz.tbi}" \
--source "reference/Homo_sapiens_assembly38.{fasta,fasta.fai}" \
--destination "reference/human_g1k_v37.{fasta,fasta.fai}" \
--chain reference/hg38ToHg19.over.chain.gz

nextflow run txt.nf -resume \
-work-dir /exports/eddie/scratch/${USER}/ad/work \
-c eddie.config \
--sumstats "liftover/*.{vcf.gz,vcf.gz.tbi}"