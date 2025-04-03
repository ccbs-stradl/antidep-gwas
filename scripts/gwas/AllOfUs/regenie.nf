/**
curl -O https://s3.amazonaws.com/plink2-assets/alpha5/plink2_linux_avx2_20240625.zip
unzip plink2_linux_avx2_20240625.zip
mkdir bin
mv plink2 bin/

curl -L -O https://github.com/rgcgithub/regenie/releases/download/v3.5/regenie_v3.5.gz_x86_64_Linux_mkl.zip
unzip regenie_v3.5.gz_x86_64_Linux_mkl.zip
mv regenie_v3.5.gz_x86_64_Linux_mkl bin/regenie

curl -s https://get.nextflow.io | bash

./nextflow run regenie.nf -resume \
-c ~/.nextflow/config -profile gls \
--bt "atc_antidep.pheno" \
--covar "cluster.covar"

**/

// analysis files
params.bt = null // binary phenotypes file
params.qt = null // quantitative phenotypes file 

params.covar = 
params.covar_list = "PC1_AVG,PC2_AVG,PC3_AVG,PC4_AVG,PC5_AVG,PC6_AVG,PC7_AVG"
params.covar_cat_list = "gender"

// genetics files
params.array = "gs://fc-aou-datasets-controlled/v7/microarray/plink_v7.1/arrays.{bed,bim,fam}"
params.srgws = "gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/acaf_threshold_v7.1/bgen/acaf_threshold.chr*.{bgen,bgen.bgi,sample}"

// qc files
params.flagged = "gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/aux/qc/flagged_samples.tsv"
params.ancestry = "gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv"

workflow {
  // ancestry assignments
  ANCESTRY_CH = Channel.fromPath(params.ancestry)

  // flagged samples
  FLAGGED_CH = Channel.fromPath(params.flagged)
  
  // parse cluster file to get names of each cluster
    CLUSTER_NAMES_CH = ANCESTRY_CH
  	.splitCsv(sep: "\t", skip: 1, header: ['iid', 'cluster', 'probs', 'features', 'cluster_other'])
  	.map { it -> it.cluster }
  	.unique()

  // array genotypes
  // assume array genotypes are FID=0, IID
  ARRAY_CH = Channel.fromFilePairs(params.array, size: 3)

	// covariates channels
	COVAR_CH = Channel
		.fromPath(params.covar)
	COVAR_LIST_CH = Channel
		.of(params.covar_list)
	COVAR_CAT_LIST_CH = Channel
	.of(params.covar_cat_list)

	/*
		Phenotypes channels
		Key input phenotype files with 'bt' or 'qt'
	*/
	
	// binary
	if(params.bt != null) {
		BT_CH = Channel
			.fromPath(params.bt, checkIfExists: true)
			.map { it -> ["bt", it] }
	} else {
		BT_CH = Channel.empty()
	}
	
	// quantitative
	if(params.qt != null) {
		QT_CH = Channel
			.fromPath(params.qt, checkIfExists: true)
			.map { it -> ["qt", it] }
	} else {
		QT_CH = Channel.empty()
	}

	// concatenate phenotype file channels
	PHENO_CH = BT_CH.mix(QT_CH)
        .map { it -> [it[0], it[1].baseName, it[1]] }	

	// set up flags required for binary or quantitative traits
	// keyed to "bt" or "qt"
	STEP1_FLAGS_CH = Channel
		.of(["bt", "--bt --minCaseCount ${params.min_cases}"],
			["qt", ""]
		)
	
	STEP2_FLAGS_CH = Channel
		.of(["bt", "--bt --af-cc --firth --approx --pThresh 0.01 --minCaseCount ${params.min_cases}"],
			["qt", ""]
		)

	FLAGS_CH = STEP1_FLAGS_CH
		.join(STEP2_FLAGS_CH)
		.map { it -> [ it[0], ["step1": it[1], "step2": it[2]] ] }


	PHENO_FLAGS_CH = PHENO_CH
		.combine(FLAGS_CH, by: 0)
		.combine(ARRAY_CH)
  	.combine(FLAGGED_CH)
  	.combine(ANCESTRY_CH)
		.combine(COVAR_CH)
		.combine(COVAR_LIST_CH)
		.combine(COVAR_CAT_LIST_CH)

	STEP1_CH = STEP1(PHENO_FLAGS_CH, CLUSTER_NAMES_CH)

	/**
		Single variants genotypes
	**/
	
	// group filenames after removing .bgen, .bgen.bgi, .sample extensions
	SRGWS_CH = Channel.fromFilePairs(params.srgws, size: 3) {it -> it.getName().minus(~/\.(bgen\.bgi|bgen|sample)/) } 

	SRGWS_PRED_CH = SRGWS_CH
		.combine(STEP1_CH)

	SRGWS_PRED_CH = SRGWS_CH
		.combine(STEP1_CH)

	STEP2_CH = STEP2(SRGWS_PRED_CH)
		.multiMap { it ->
			regenie: it[1]
			log: it[2]
		}

	// // group per-chromosome sumstats for each phenotype*cluster
	STEP2_REGENIE_CH = STEP2_CH.regenie
		.flatten()
	  	.map { it -> [(it.name =~ /step2_(.+)\.chr[0-9XY]+_(.+)\.regenie\.gz/)[0][1..2], it] }	
		.groupTuple()
		.map { it -> [it[0][0], it[0][1], it[1]] }

	GWAS_CH = MERGE(STEP2_REGENIE_CH)

}

/** genome regression
	  convert input files to FID=0, IID
		QC participants and SNPs for each cluster
		run genome regression
**/
process STEP1 {
	tag "${pheno}-${cluster}"

	machineType 'c2d-highcpu-8'
  container "${System.getenv('ARTIFACT_REGISTRY_DOCKER_REPO')}/r-base:4.4.1"
	//scratch true, mode: 'copy'

	input:
	tuple val(traits), val(pheno), path(phenos), val(flags), val(array), path(bedbimfam), path(flagged), path(ancestry), path(covars), val(covar_list), val(covar_cat_list)
	each cluster

	output:
	tuple val(traits), val(pheno), path(phenos, includeInputs: true), path(covars, includeInputs: true), val(covar_list), val(covar_cat_list), val(cluster), val(flags), val("step1_${traits}-${pheno}-${cluster}"), path("step1_*") 


	script:
	"""
	cat ${phenos} | awk -v OFS='\t' '{if(NR == 1){print \$0} else {\$1 = 0; print \$0}}' > ${phenos.baseName}.fid0.pheno
	cat ${covars} | awk -v OFS='\t' '{if(NR == 1){print \$0} else {\$1 = 0; print \$0}}' > ${covars.baseName}.fid0.covar

  cat ${flagged} | awk -v OFS='\t' '{if(NR == 1) {print "FID", "IID"} else {print 0, \$1}}' > remove.ids
  cat ${ancestry} | awk -v OFS='\t' '{if(NR == 1) {print "FID", "IID", "cluster"} else {print 0, \$1, \$2}}' > clusters.txt

    plink2 \
      --bfile ${array} \
      --maf 0.01 --mac 100 --geno 0.1 --hwe 1e-15 \
      --thin-count 500000 \
      --mind 0.1 \
      --keep-col-match clusters.txt ${cluster} \
      --keep-col-match-num 3 \
      --remove remove.ids \
      --write-snplist \
      --write-samples --no-id-header \
      --out ${cluster} \
      --threads 8 \
      --memory 16000

	regenie --step 1 ${flags.step1} \
	--phenoFile ${phenos.baseName}.fid0.pheno \
	--bed ${array} \
	--extract ${cluster}.snplist \
	--keep ${cluster}.id \
	--covarFile ${covars.baseName}.fid0.covar \
	--covarColList ${covar_list} \
	--catCovarList ${covar_cat_list} \
	--out step1_${traits}-${pheno}-${cluster} \
	--use-relative-path \
	--bsize 1000 \
	--threads 8
	"""
}

/** genome association
	  convert input files to FID=0, IID
		run genome association
**/
process STEP2 {
	tag "${pheno}-${cluster}-${srgws}"

	machineType 'c2d-highcpu-8'
  container "${System.getenv('ARTIFACT_REGISTRY_DOCKER_REPO')}/r-base:4.4.1"
	//scratch true, mode: 'copy'
	errorStrategy 'finish'

	input:
	tuple val(srgws), path(bgen), val(traits), val(pheno), path(phenos), path(covars), val(covar_list), val(covar_cat_list), val(cluster), val(flags), val(step1), path(pred_list) 

	output:
	tuple val("step2_${pheno}-${cluster}-${srgws}"), path("*.regenie.gz"), path("*.log")

	script:
	"""
	cat ${phenos} | awk -v OFS='\t' '{if(NR == 1){print \$0} else {\$1 = 0; print \$0}}' > ${phenos.baseName}.fid0.pheno
	cat ${covars} | awk -v OFS='\t' '{if(NR == 1){print \$0} else {\$1 = 0; print \$0}}' > ${covars.baseName}.fid0.covar
	cat ${srgws}.sample | awk -v OFS='\t' '{if(NR == 1){print \$0} else {\$1 = 0; print \$0}}' > ${srgws}.fid0.sample

	regenie \
	--step 2 ${flags.step2} \
	--bgen ${srgws}.bgen \
	--sample ${srgws}.fid0.sample \
	--phenoFile ${phenos.baseName}.fid0.pheno \
	--covarFile ${covars.baseName}.fid0.covar \
	--covarColList ${covar_list} \
	--catCovarList ${covar_cat_list} \
	--pred ${step1}_pred.list \
	--bsize 400 \
	--minMAC 40 \
	--minINFO 0.1 \
	--out step2_${pheno}-${cluster}-${srgws} \
	--gz \
	--threads 8 
"""
}


/** merge sumstats 
		order by chromosome and position
		name chromosomes as chrM
**/
process MERGE {
	tag "${dataset}-${pheno}"
	publishDir "sumstats/${dataset}", mode: 'copy'
        	
	machineType 'c4-highmem-4'
  container "${System.getenv('ARTIFACT_REGISTRY_DOCKER_REPO')}/rocker/tidyverse"
	errorStrategy 'finish'
	input:
	tuple val(dataset), val(pheno), path(regenies)
	
	output:
	tuple val(dataset), val(pheno), path("*.regenie.gz")
	
	script:
	"""
	#!Rscript
	library(dplyr)
	library(readr)
	library(stringr)
	
	regenie_paths <- str_split("${regenies}", pattern = ' ')[[1]]
	
	regenies <-
	bind_rows(
		lapply(regenie_paths,
			   read_table,
			   col_types=cols(CHROM=col_character(),
				       EXTRA=col_character()),
			   na = c('', 'NA', '-nan' )
		)
	) |>
	mutate(chr = case_match(CHROM,
			"X" ~ 23,
			"Y" ~ 24,
			.default = as.integer(CHROM))) |>
	arrange(chr, GENPOS) |>
	mutate(CHROM = str_c("chr", CHROM)) |>
	select(-chr)
	
	write_tsv(regenies, "${dataset}-${pheno}.regenie.gz")
	"""
}
