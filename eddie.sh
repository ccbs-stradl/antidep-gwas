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
nextflow run workflows/format_gwas.nf -resume \
-work-dir /exports/eddie/scratch/${USER}/ad/work/format \
-c config/eddie.config

# convert hg38 sumstats to VCF
nextflow run workflows/vcf.nf -resume \
--sumstats "results/format/gwas/GRCh38/*.{txt,json,csv}" \
--chr 'chr#' \
--publish "results/vcf/gwas/GRCh38" \
--assembly "reference/Homo_sapiens_assembly38.{fasta,fasta.fai,dict}" \
--dbsnp "reference/dbsnp.v153.hg38.vcf.{gz,gz.tbi}" \
-work-dir /exports/eddie/scratch/${USER}/ad/work/vcf \
-c config/eddie.config

# convert hg19 sumstats to VCF
nextflow run workflows/vcf.nf -resume \
--sumstats "results/format/gwas/GRCh37/*.{txt,json,csv}" \
--chr '#' \
--publish "results/vcf/gwas/GRCh37" \
--assembly "reference/human_g1k_v37.{fasta,fasta.fai,dict}" \
--dbsnp "reference/dbsnp.v153.b37.vcf.{gz,gz.tbi}" \
-work-dir /exports/eddie/scratch/${USER}/ad/work/vcf \
-c config/eddie.config

# liftover hg19 VCFs to hg38
nextflow run workflows/liftover.nf -resume \
--sumstats "results/vcf/gwas/GRCh37/*.{vcf.gz,vcf.gz.tbi}" \
--source "reference/human_g1k_v37.{fasta,fasta.fai}" \
--destination "reference/Homo_sapiens_assembly38.{fasta,fasta.fai}" \
--chain reference/hg19ToHg38.over.chain.gz \
--publish results/liftover/gwas \
-work-dir /exports/eddie/scratch/${USER}/ad/work/liftover \
-c config/eddie.config 

# move lifted over file
mv results/liftover/gwas/*.Homo_sapiens_assembly38.vcf.* results/vcf/gwas/GRCh38
rename Homo_sapiens_assembly38.vcf vcf results/vcf/gwas/GRCh38/*.Homo_sapiens_assembly38.vcf.*
cp results/vcf/gwas/GRCh37/*.csv results/vcf/gwas/GRCh38/

###
### Run meta analysis
###

# run meta-analysis
nextflow run workflows/meta.nf -resume \
-work-dir /exports/eddie/scratch/${USER}/ad/work/meta \
-c config/eddie.config

# convert meta-analysis to VCF
nextflow run workflows/format_meta.nf -resume \
-work-dir /exports/eddie/scratch/${USER}/ad/work/meta \
-c config/eddie.config

nextflow run workflows/vcf.nf -resume \
--sumstats "format/meta/GRCh38/*.{txt,json,csv}" \
--chr 'chr#' \
--publish "vcf/meta/GRCh38" \
--assembly "reference/Homo_sapiens_assembly38.{fasta,fasta.fai,dict}" \
--dbsnp "reference/dbsnp.v153.hg38.vcf.{gz,gz.tbi}" \
-work-dir /exports/eddie/scratch/${USER}/ad/work/vcf \
-c config/eddie.config -with-trace

# liftover hg19 VCFs to hg38
nextflow run workflows/liftover.nf -resume \
--sumstats "vcf/meta/GRCh38/*.{vcf.gz,vcf.gz.tbi}" \
--source "reference/Homo_sapiens_assembly38.{fasta,fasta.fai}" \
--destination "reference/human_g1k_v37.{fasta,fasta.fai}" \
--chain reference/hg38ToHg19.over.chain.gz \
--publish liftover/meta \
-work-dir /exports/eddie/scratch/${USER}/ad/work/liftover \
-c config/eddie.config 

# move lifted over file and copy sidecar files
mv liftover/meta/*-fixed-*.human_g1k_v37.vcf.gz.* vcf/meta/GRCh37
rename human_g1k_v37.vcf vcf vcf/meta/GRCh38/*.human_g1k_v37.vcf.*
cp vcf/meta/GRCh37/*.csv vcf/meta/GRCh38/

###
### Run downstream analyses
###

# Run mBAT-combo on build hg19/GRCh37
nextflow run workflows/genes.nf -resume \
-work-dir /exports/eddie/scratch/${USER}/ad/work_hg19 \
-c config/eddie.config \
--build 'hg19'

# Run mBAT-combo on build hg38/GRCh38
nextflow run workflows/genes.nf -resume \
-work-dir /exports/eddie/scratch/${USER}/ad/work_hg38 \
-c config/eddie.config \
--build 'hg38'

# Run popcorn
nextflow run workflows/popcorn.nf -resume \
-work-dir /exports/eddie/scratch/${USER}/ad/work_hg38 \
-c config/eddie.config

nextflow run workflows/txt.nf -resume \
-work-dir /exports/eddie/scratch/${USER}/ad/work \
-c config/eddie.config \
--sumstats "liftover/*.{vcf.gz,vcf.gz.tbi}"

# Run SuSiEx on build hg19
nextflow run workflows/fine_mapping.nf -resume \
-work-dir /exports/eddie/scratch/${USER}/ad/work \
-c config/eddie.config -with-dag fineMapping/fine_mapping_dag.png

# Plot results of SuSiEx
Rscript fine_mapping_plots.R

#
