/* GWAS in regenie with TOPMed imputation (GRCh38) */

nextflow.enable.dsl=2

params.bt = null // binary phenotypes file
params.qt = null // quantitative phenotypes file 
params.keep = "genscot/genetics/genotypes/GS20K_PLINK_files/QCd_data/QCdGS20K.keep"
params.remove = "genscot/phenotypes/overlaps/GSindividualsOverlappingWithUKB500K.txt"
params.bfile = "genscot/genetics/genotypes/GS20K_PLINK_files/QCd_data/QCdGS20K.{bed,bim,fam}"
params.bgen = "TOPMed/TOPMedFreeze5b/GS20K/nomono_I1/bgen/GS20K_chr*_TOPMedFreeze5b_nomono_I1.{bgen,bgen.bgi,sample}"
params.covar = "genscot/phenotypes/agesex_pcs1-4.txt"
params.covar_list = "age,C1,C2,C3,C4"
params.covar_cat_list = "sex"

// reference file parameters
params.dbsnp = "resources/dbSNP/dbsnp.v153.hg38.vcf.{gz,gz.tbi}"


workflow {

/*
	Inputs
*/
	// participants to analyze
	KEEP_CH = Channel
		.fromPath(params.keep, checkIfExists: true)
		
	// participants to not analyze
	REMOVE_CH = Channel
		.fromPath(params.remove, checkIfExists: true)
		
	// plink genotypes
	BFILE_CH = Channel
		.fromFilePairs(params.bfile, size: 3, checkIfExists: true)
		
	// bgen imputed
	BGEN_CH = Channel
		.fromFilePairs(params.bgen, size: 3, checkIfExists: true) { it -> it.simpleName }
		
	// dbSNP sites references file
	DBSNP_CH = Channel
		.fromFilePairs(params.dbsnp, size: 2)
		
	// covariates channels
	COVAR_CH = Channel
		.fromPath(params.covar, checkIfExists: true)
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
	
/*
		Preprocessing
*/
		
	// genotyped SNPs
	SNP_CH = BFILE_SNPS(BFILE_CH)
		
	// select SNPs from the references file
	DBSNP_SELECT_CH = SELECT(DBSNP_CH, SNP_CH)
	
	// update genotype map to affect liftover, QC samples and variants
	BFILE_OVER_CH = UPDATE_MAP(BFILE_CH, DBSNP_SELECT_CH, KEEP_CH, REMOVE_CH)
		
	// set up flags required for binary or quantitative traits
	// keyed to "bt" or "qt"
	STEP1_FLAGS_CH = Channel
		.of(["bt", "--bt --minCaseCount 200"],
			["qt", ""]
		)
	
	STEP2_FLAGS_CH = Channel
		.of(["bt", "--bt --af-cc --firth --approx --pThresh 0.01"],
			["qt", ""]
		)
	
	// join phenotype files with flags using bt/qt key,
	// then tack on all other support files
	PHENO_FLAGS_BFILE_CH = PHENO_CH
		.join(STEP1_FLAGS_CH)
		.combine(BFILE_OVER_CH)
		.combine(COVAR_CH)
		.combine(COVAR_LIST_CH)
		.combine(COVAR_CAT_LIST_CH)
	
/*
		GWAS
*/
		
	// step 1 loco for each phenotype file
	// returns pred.list keyed to bt/qt
	STEP1_CH = STEP1(PHENO_FLAGS_BFILE_CH)
	
	// using bt/qt key to join phenotype files with
	// step2 flags and pred.list from step 1
	STEP1_FLAGS_PFILE_CH = PHENO_CH
		.join(STEP2_FLAGS_CH)
		.join(STEP1_CH)
		.combine(BGEN_CH)
		.combine(COVAR_CH)
		.combine(COVAR_LIST_CH)
		.combine(COVAR_CAT_LIST_CH)
	
	// variant association step, then split apart .regenie and .log files
	STEP2_CH = STEP2(STEP1_FLAGS_PFILE_CH)
		.multiMap { it ->
			regenie: it[0]
			log: it[1]
		}
		
/*
		Merging and output
*/
	// group per-chromosome sumstats for each phenotype
	STEP2_REGENIE_CH = STEP2_CH
		.regenie
		.flatten()
		
	STEP2_LOG_CH =
	STEP2_CH
		.log
		.flatten()
		.map { it -> [(it.name =~ /step2_([bq]t)_(.+)_chr[0-9X]+_(.+)--\.log/)[0][1..3], it]}
		.groupTuple()
		
	// get names of keep, and remove files to make name of sumstats file
	KEEP_VAL =
	KEEP_CH
		.map { it -> it.baseName }
		.first()
		
	REMOVE_VAL =
	REMOVE_CH
		.map { it -> it.baseName }
		.first()

	// get rsIDs from sites reference
	// tack on topmed ref
	GWAS_REF_CH = STEP2_REGENIE_CH
		.combine(DBSNP_CH)
	
	REGIONS_CH = REGIONS(GWAS_REF_CH)
	
	GWAS_RSID_CH = CPID_RSID(REGIONS_CH)
		.flatten()
		.map { it -> [(it.name =~ /step2_([bq]t)_(.+)_chr[0-9X]+_(.+)--_(.+)\.regenie/)[0][1..4], it]}
		.groupTuple()
	
	// merge chromosomes together for each GWAS, make final filename
	GWAS_CH = MERGE(GWAS_RSID_CH, KEEP_VAL, REMOVE_VAL)

	// merge logs
	LOGS_CH = LOGS(STEP2_LOG_CH, KEEP_VAL, REMOVE_VAL)

	// make plots
	MANH_CH = MANHATTAN(GWAS_CH) 
	
}

/* Get SNP list in the genotype bim file */
process BFILE_SNPS {
	tag "${dataset}"
	
	executor 'local'
	
	cpus = 1
	memory = 1.GB
	time = '1m'
	
	input:
	tuple val(dataset), path(bfile)
	
	output:
	path("${dataset}.list")
	
	script:
	"""
	cat ${dataset}.bim | awk '{print \$2}' > ${dataset}.list
	"""
}

/* Select a SNP list from the dbSNP VCF sites file */
process SELECT {
	tag "${dbsnp}"
	
	cpus = 1
	memory = 8.GB
	time = '30m'
	
	input:
	tuple val(dbsnp), path(vcf)
	path(list)
	
	output:
	path("${vcf[0].simpleName}-select.vcf.gz")
	
	script:
	"""
	gatk SelectVariants \
	-V ${dbsnp}.gz \
	--keep-ids ${list} \
	-O ${vcf[0].simpleName}-select.vcf.gz
	"""
}

/* Liftover matched SNPs by updating the bfile map 
   Keep and removing specified samples
   Update IDs to be double IIDs to match BGEN files */
process UPDATE_MAP {
	tag "${cohort}"
	
	cpus = 1
	memory = 8.GB
	time = '30m'
	
	input:
	tuple val(cohort), path(bfile)
	path(vcf)
	path(keep)
	path(remove)
	
	output:
	tuple val("${cohort}-lifted"), path("${cohort}-lifted.{bed,bim,fam}")
	
	script:
	"""
	gunzip -c ${vcf} | grep -v '^#' | awk '{print \$3}' > snps.list
	gunzip -c ${vcf} | grep -v '^#' | awk '{print \$3, \$2}' > snps.map
	cat ${cohort}.fam | awk '{print \$1, \$2, \$2, \$2}' > double.ids
	
	plink2 --bfile ${cohort} \
	--make-pgen \
	--keep ${keep} \
	--remove ${remove} \
	--extract snps.list \
	--update-map snps.map \
	--sort-vars \
	--out ${cohort}-lifted \
	--threads ${task.cpus} \
	--memory ${task.memory.bytes.intdiv(1000000)}
	
	plink2 --pfile ${cohort}-lifted \
	--make-bed \
	--update-ids double.ids \
	--maf 0.01 --mac 100 --geno 0.1 --hwe 1e-15 \
	--out ${cohort}-lifted \
	--threads ${task.cpus} \
	--memory ${task.memory.bytes.intdiv(1000000)}
	"""
}

/* regenie step1 loco */
process STEP1 {
	tag "${traits}"

	cpus = 8
	memory = 16.GB
	time = '1h'

	input:
	tuple val(traits), path(pheno), val(flags), val(bfile_prefix), path(bfile), path(covar), val(covar_list), val(covar_cat_list)

	output:
	tuple val(traits), path("step1_${traits}*")

	script:
	"""
	doubleid.R ${pheno} iid-${pheno}
	doubleid.R ${covar} iid-${covar}
	
	regenie --step 1 ${flags} \
	--phenoFile iid-${pheno} \
	--bed ${bfile_prefix} \
	--covarFile iid-${covar} \
	--covarColList ${covar_list} \
	--catCovarList ${covar_cat_list} \
	--out step1_${traits} \
	--bsize 1000 \
	--threads ${task.cpus}
	"""

}


/* regenie step 2 association test */
process STEP2 {
	tag "Step 2 ${traits}"
	
	cpus = 8
	memory = 16.GB
	time = '24h'
	
	input:
	tuple val(traits), path(pheno), val(flags), path(step1), val(bgen_prefix), path(bgen), path(covar), val(covar_list), val(covar_cat_list)
	
	output:
	tuple path("step2_${traits}_${bgen_prefix}*.regenie"), path("step2_${traits}_${bgen_prefix}*.log")
	
	script:
	"""
	doubleid.R ${pheno} iid-${pheno}
	doubleid.R ${covar} iid-${covar}
	
	regenie \
	  --step 2 ${flags} \
	  --bgen ${bgen_prefix}.bgen \
	  --sample ${bgen_prefix}.sample \
	  --phenoFile iid-${pheno} \
	  --covarFile iid-${covar} \
	  --covarColList ${covar_list} \
	  --catCovarList ${covar_cat_list} \
	  --pred step1_${traits}_pred.list \
	  --bsize 400 \
	  --minMAC 20 \
	  --minINFO 0.1 \
	  --out step2_${traits}_${bgen_prefix}-- \
	  --threads ${task.cpus}
	"""
	

}


/* get RSIDs matching the all GWAS CPIDs from the sites reference */
process REGIONS {
	tag "${regenie.simpleName}"
	
	cpus = 1
	memory = 4.GB
	time = '3h'
	
	input:
	tuple path(regenie), val(dbsnp), path(dbsnp_vcf)
	
	output:
	tuple path(regenie, includeInputs: true), path("${dbsnp_vcf[0].simpleName}-rsids.tsv.gz")
	
	script:
	"""
	# get chrom and position in 0-indexed bed format
	cat ${regenie} | awk '{OFS="\\t"; if(NR > 1) {sub("23", "X", \$1); print "chr"\$1, \$2-1, \$2}}' > targets.bed
	# block zip and index
	bgzip targets.bed
	tabix targets.bed.gz

	bcftools view --regions-file targets.bed.gz --targets-file targets.bed.gz ${dbsnp}.gz |\
	bcftools query -f "%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\n" |\
	split_alts.awk |\
	gzip -c > ${dbsnp_vcf[0].simpleName}-rsids.tsv.gz
	"""
}

/* rename CPID to rsID */
process CPID_RSID {
	tag "${regenie}"
		
	cpus = 1
	memory = 16.GB
	time = '1h'
	
	input:
	tuple path(regenie), path(dbsnp)
	
	output:
	path "*.tsv.gz"
	
	script:
	"""
	#!Rscript
	library(dplyr)
	library(readr)
	library(stringr)
	library(fastman)
	
	regenie <- read_table("${regenie}", col_types = cols(GENPOS=col_integer())) 
	dbsnp <- read_tsv("${dbsnp}", comment = "##",
		col_names=c("CHROM", "GENPOS", "RSID", "REF", "ALT"),
		col_types = cols(CHROM = col_character(), GENPOS=col_integer()))
		
	dbsnp_count <- dbsnp |>
		group_by(CHROM, GENPOS, REF, ALT) |>
		mutate(n = n()) |>
		ungroup()
		
	dbsnp_dedup <- dbsnp_count |>
		filter(n == 1) |>
		select(-n)
		
	regenie_rsids <- regenie |>
		mutate(CHR = CHROM) |>
		mutate(CHROM = str_c("chr", if_else(CHR == 23, true = "X", false = as.character(CHR)))) |>
		mutate(CPID = if_else(ID == ".", true=str_c(CHROM, ":", GENPOS, ":", ALLELE0, ":", ALLELE1), false=ID)) |>
		left_join(dbsnp_dedup, by = c("CHROM", "GENPOS", "ALLELE0"="ALT", "ALLELE1"="REF")) |>
		mutate(ID = coalesce(RSID, CPID)) |>
		select(-RSID, -CPID)

	write_tsv(regenie_rsids, "${regenie}.tsv.gz")
	"""
}

/* merge sumstats */
process MERGE {
	tag "merging ${pheno[3]}"
	
	publishDir "sumstats", mode: 'copy'
	
	cpus = 1
	memory = 16.GB
	time = '10m'
	
	input:
	tuple val(pheno), path(regenies)
	val keep
	val remove
	
	output:
	path "*.regenie.gz"
	
	script:
	"""
	#!Rscript
	library(dplyr)
	library(readr)
	library(stringr)
	library(glue)
	
	regenie_paths <- str_split("${regenies}", pattern = ' ')[[1]]

	regenies <-
	bind_rows(
		lapply(regenie_paths, read_table)
	) |>
	arrange(CHR, GENPOS) |>
	select(-CHR)
	
	write_tsv(regenies, "${pheno[1]}_${pheno[2]}-${pheno[0]}-${pheno[3]}-${keep}-no${remove}.regenie.gz")
	"""
}

/* merge logs */
process LOGS {
	tag "logs ${traits[0]}"
	publishDir 'sumstats', mode: 'copy'
	
	executor 'local'
	
	cpus = 1
	memory = 16.GB
	time = '10m'
	
	input:
	tuple val(traits), path(logs)
	val keep
	val remove
	
	output:
	path "*.log"
	
	script:
	"""
	cat ${logs} > "${traits[1]}-${traits[2]}-${traits[0]}-${keep}-no${remove}.log"
	"""
}



/* manhattan plot */
process MANHATTAN {
	tag "${regenie}"
	publishDir 'sumstats', mode: 'copy'
	
	cpus = 1
	memory = 16.GB
	time = '20m'
	
	input:
	path(regenie)
	
	output:
	path "*.png"
	
	script:
	"""
	#!Rscript
	library(dplyr)
	library(readr)
	library(fastman)
	library(stringr)
	
	regenie <- read_tsv("${regenie}", col_types=cols(CHROM = col_character())) |>
		mutate(P=10^(-LOG10P)) |>
		mutate(CHROM = if_else(CHROM == "chrX", true = 23L, false = as.integer(str_replace(CHROM, pattern="chr", replacement=""))))
	
	png("${regenie.baseName}.manhattan.png", width=10, height=6, units="in", res=300)
	fastman(regenie, chr = "CHROM", bp = "GENPOS", p = "P", maxP=NULL)
	dev.off()
	
	png("${regenie.baseName}.qq.png", width=6, height=6, units="in", res=300)
	fastqq(p1=pull(regenie, P), maxP=NULL)
	dev.off()
	"""
}