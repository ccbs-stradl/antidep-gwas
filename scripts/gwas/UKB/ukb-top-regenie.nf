/* UKB GWAS in regenie on RAP with TOPMed imputation (GRCh38) */


nextflow.enable.dsl=2

// nextflow run ukb-top-regenie.nf -resume -config eddie.config -work-dir /exports/eddie/scratch/$USER/ukb/top

params.bt = null // binary phenotypes file
params.qt = null // quantitative phenotypes file 
params.keep = "/exports/igmm/eddie/GenScotDepression/data/ukb/genetics/input_filters/similarity_clusters/rf_hgdp1kg_clusters.keep"
params.remove = "/exports/igmm/eddie/GenScotDepression/data/ukb/genetics/input_filters/v2/PGC.remove"
params.bfile = "/exports/igmm/eddie/GenScotDepression/data/ukb/genetics/impv2/bfile/autosome-grch38/ukb_v2-grch38.{bed,bim,fam}"
params.clusters = "/exports/igmm/eddie/GenScotDepression/data/ukb/genetics/input_filters/similarity_clusters/ukb_randomforest_clusters.tsv"
params.covar = "/exports/igmm/eddie/GenScotDepression/data/ukb/genetics/input_filters/similarity_clusters/ukb_randomforest_clusters.covar"
params.covar_list = "PC1,PC2,PC3,PC4,PC5,PC6"
params.covar_cat_list = "sex,genotyping"
params.min_cases = 80

// reference file parameters
params.dbsnp = "/exports/igmm/eddie/GenScotDepression/data/resources/dbSNP/dbsnp.v153.hg38.vcf.{gz,gz.tbi}"

workflow {
	// participants to analyze
	KEEP_CH = Channel
		.fromPath(params.keep, checkIfExists: true)
		
	// participants to not analyze
	REMOVE_CH = Channel
		.fromPath(params.remove, checkIfExists: true)
		
	// plink genotypes
	BFILE_CH = Channel
		.fromFilePairs(params.bfile, size: 3, checkIfExists: true)
		
	// genetic similarity clusters
	CLUSTERS_CH = Channel
		.fromPath(params.clusters, checkIfExists: true)
		
	// parse cluster file to get names of each cluster
	CLUSTER_NAMES_CH = CLUSTERS_CH
		.splitCsv(sep: "\t", skip: 1, header: ['fid', 'iid', 'cluster'])
		.map { it -> it.cluster }
		.unique()

	// covariates channels
	COVAR_CH = Channel
		.fromPath(params.covar, checkIfExists: true)
	COVAR_LIST_CH = Channel
		.of(params.covar_list)
	COVAR_CAT_LIST_CH = Channel
	.of(params.covar_cat_list)

	// dbSNP  references file
	DBSNP_CH = Channel
		.fromFilePairs(params.dbsnp, size: 2)
	
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
		

	/*
		Analysis processes
		Key files with trait ('bt' or 'qt') and cluster
	*/

	/*
		Local compute
	*/

	// Genotype QC
	BFILE_CLUSTERS_CH = BFILE_CH
		.combine(KEEP_CH)
		.combine(REMOVE_CH)
		.combine(CLUSTERS_CH)
		
	QC_CH = QC(BFILE_CLUSTERS_CH, CLUSTER_NAMES_CH)
	
	// join phenotype files with flags using bt/qt key,
	// then tack on all other support files
	PHENO_FLAGS_BFILE_CH = PHENO_CH
		.combine(FLAGS_CH, by: 0)
		.combine(BFILE_CH)
		.combine(QC_CH)
		.combine(COVAR_CH)
		.combine(COVAR_LIST_CH)
		.combine(COVAR_CAT_LIST_CH)

	// step 1 loco for each phenotype file
	// returns pred.list keyed to bt/qt
	STEP1_CH = STEP1(PHENO_FLAGS_BFILE_CH)
		.combine(COVAR_CH)
		.combine(COVAR_LIST_CH)
		.combine(COVAR_CAT_LIST_CH)

	/*
		RAP compute
	*/

    // upload files to dnanexus
    STAGE_CH = STAGE(STEP1_CH)

    // get files list for input into step 2
    IIN_CH = LS(STAGE_CH)
        .map { it -> [it[0], it[1], it[2].splitCsv()] }
        .transpose()
        .map { it -> [it[0], "-iin=\"/work/${it[0]}/${it[1]}/${it[2][0]}\""]}
        .groupTuple()

    STAGE_IIN_CH = STAGE_CH
        .join(IIN_CH)

    // // Run step 2 on each chromosome
	CHR_CH = Channel.of(1..22, 'X')

	// Check job status
	FIND_CH = FIND(STAGE_IIN_CH, CHR_CH)
		.map { it -> [it, file(it[3]).readLines()]}
	
	// Branch depending on if job has been submitted
	JOBS_CH = FIND_CH
		.branch {
			todo: it[1].empty
			monitor: !it[1].empty
		}

	// Jobs to be submitted to dx run
	RUN_CH =
	JOBS_CH
		.todo
		.map { it -> it[0] }

	// Run Step 2 on dx
    STEP2_CH = STEP2(RUN_CH)

	// Running and submitted jobs to monitor
	MONITOR_CH =
	JOBS_CH
		.monitor
		.map { it -> [it[0][0], it[0][1], it[0][2], it[1][0]]}
		.concat(STEP2_CH)


	// Get job meta data from JSON
	def slurp = new groovy.json.JsonSlurper()

	DESCRIBE_CH = DESCRIBE(MONITOR_CH)
		.map { it -> [it[0], it[1], it[2], it[3], slurp.parseText(file(it[4]).text)]}

	/*
	def now = new Date()
	DESCRIBE_CH
		.map { it -> [ "state": it[4]['state'], 'dataset': it[0], 'hash': it[1], 'chr': it[2], 'jobid': it[3], 'start': it[4]['startedRunning'].toBigInteger(), 'end': it[4]['stoppedRunning'] ? it[4]['stoppedRunning'].toBigInteger() : now.getTime() ] }
		.view { it -> "${it['state']}: ${it['dataset']}, ${it['chr']}, ${it['jobid']}, ${(it['end'] - it['start']).millis}"}
	*/

	// // Clean up failed jobs and get files for finished jobs
	STATUS_CH = DESCRIBE_CH
		.branch {
			failed: it[4]['state'] == 'failed'
			done: it[4]['state'] == 'done'
		}

	// clean up metadata properties from failed jobs
	FAILED_CH = STATUS_CH.failed
	CLEAN(FAILED_CH)

    // download regenie files
	DONE_CH = STATUS_CH.done
		.map { it -> [[it[0], it[1]], ["state": it[4]['state'], 'dataset': it[0], 'hash': it[1], 'chr': it[2], 'jobid': it[3], 'start': it[4]['startedRunning'], 'end': it[4]['stoppedRunning'], 'price': it[4]['totalPrice'] ] ] }
	
	DATASET_DONE_CH = DONE_CH
		.groupTuple()
		.map { it -> [it[0], it[1].size() ] }
		.filter { it -> it[1] == 23 }
		.map { it -> [it[0][0], it[0][1]]}

    FETCH_CH = FETCH(DATASET_DONE_CH)
		.multiMap { it ->
			regenie: it[1]
			log: it[2]
		}
	
	/*
		Local compute
	*/

	REGENIE_CH = FETCH_CH
		.regenie
		.flatten()
		.map { it -> [(it.name =~ /step2_(.+)_chr.+_(.+)\.regenie/)[0][1..2], it]}
		.groupTuple()
		.map { it -> [it[0][0], it[0][1], it[1]] }

	KEEP_VAL =
	KEEP_CH
		.map { it -> it.baseName }
		.first()
		
	REMOVE_VAL =
	REMOVE_CH
		.map { it -> it.baseName }
		.first()

	GWAS_CH = MERGE(REGENIE_CH, KEEP_VAL, REMOVE_VAL)
	MANH_CH = MANHATTAN(GWAS_CH)

	// get rsIDs from TOPMed sites reference
	// collect all sumstats together, then 
	// tack on topmed ref
	GWAS_REF_CH = GWAS_CH
		.combine(DBSNP_CH)

	REGIONS_CH = REGIONS(GWAS_REF_CH)
	
	GWAS_RSID_CH = CPID_RSID(REGIONS_CH)
	

	// Summarize runtime
	RUNTIME_CH =
	DONE_CH
	  .map { it -> [it[1].dataset, it[1].hash, it[1].chr, it[1].jobid, it[1].start, it[1].end, it[1].price].join(',') }
	  .collectFile(newLine: true)

	RUNTIME(RUNTIME_CH)
	
}

/* Cluster genotypes */
process QC {
	tag "QCing genotypes ${cluster}"
	
	cpus = 1
	memory = 8.GB
	time = '30m'
	
	input:
	tuple val(bfile), path(bedbimfam), path(keep), path(remove), path(clusters)
	each cluster
	
	
	output:
	tuple val(cluster), path("${cluster}.snplist"), path("${cluster}.id")
	
	script:
	"""
	plink2 \
	  --bfile ${bfile} \
	  --maf 0.01 --mac 100 --geno 0.1 --hwe 1e-15 \
	  --mind 0.1 \
	  --keep ${keep} \
	  --keep-col-match ${clusters} ${cluster} \
	  --keep-col-match-num 3 \
	  --remove ${remove} \
	  --write-snplist \
	  --write-samples --no-id-header \
	  --out ${cluster} \
	  --threads ${task.cpus} \
	  --memory ${task.memory.bytes.intdiv(1000000)}
	"""
}

/* regenie step1 loco */
process STEP1 {
	tag "${cluster}-${traits}-${pheno.baseName}"
	
	// processes might not produce any output if there are not
	// enough cases for any phenotype
	//errorStrategy 'ignore'

	cpus = 8
	memory = 64.GB
	time = '6h'

	input:
	tuple val(traits), path(pheno), val(flags), val(bfile), path(bedbimbam), val(cluster), path(extract), path(keep), path(covar), val(covar_list), val(covar_cat_list)

	output:
	tuple val(traits), val(cluster), path(pheno, includeInputs: true), val(flags), val("step1_${cluster}-${traits}-${pheno.baseName}"), path("step1_${cluster}-${traits}-${pheno.baseName}*")

	script:
	"""
	regenie --step 1 ${flags.step1} \
	--phenoFile ${pheno} \
	--bed ${bfile} \
	--extract ${extract} \
	--keep ${keep} \
	--covarFile ${covar} \
	--covarColList ${covar_list} \
	--catCovarList ${covar_cat_list} \
	--out step1_${cluster}-${traits}-${pheno.baseName} \
    --use-relative-path \
	--bsize 1000 \
	--threads ${task.cpus}
	"""

}

/* Stage to DNAnexus 
   Get hash to keep track of these inputs
   Upload files to project
*/
process STAGE {
    tag "${cluster}-${traits}-${pheno.baseName}"

    executor 'local'

    cpus = 1
    memory = 4.GB
    time = '10m'

    input:
    tuple val(traits), val(cluster), path(pheno), val(flags), val(step1), path(loco), path(covar), val(covar_list), val(covar_cat_list)

    output:
    tuple val("${cluster}-${traits}-${pheno.baseName}"), env(hash), val(traits), val(cluster), val(pheno.name), val(covar.name), val(flags.step2), val("${step1}_pred.list"), val(covar_list), val(covar_cat_list)

    script:
    """
	hash=\$(basename \$(pwd))
    dx upload --parents --wait --no-progress --singlethread --destination="/work/${cluster}-${traits}-${pheno.baseName}/\${hash}" ${pheno} ${covar} ${loco}
    """
}

/* list files that are required in the workflow */
process LS {
    tag "${dataset}"

    executor 'local'

    cpus = 1
    memory = 4.GB
    time = '10m'

    input:
    tuple val(dataset), val(hash), val(traits), val(cluster), val(pheno), val(covar), val(flags), val(pred), val(covar_list), val(covar_cat_list)

    output:
    tuple val(dataset), val(hash), path("files.txt")

    script:
    """
    dx ls "/work/${dataset}/${hash}" > files.txt
    """

}

/* Find job status */
process FIND {
	tag "${dataset} ${chr}"

	executor 'local'
	cache true

	cpus = 1
	memory = 4.GB
	time = '10m'

	input:
	tuple val(dataset), val(hash), val(traits), val(cluster), val(pheno), val(covar), val(flags), val(pred), val(covar_list), val(covar_cat_list), val(iin)
	each chr

	output:
	tuple val(dataset), val(hash), val(chr), path("jobid.txt"), val(traits), val(cluster), val(pheno), val(covar), val(flags), val(pred), val(covar_list), val(covar_cat_list), val(iin)

	script:
	"""
	dx find jobs --brief --property hash=${hash} --property chr=${chr} > jobid.txt
	"""


}

/* Get job runtime details */
process DESCRIBE {
	tag "${dataset} ${chr}"
	publishDir 'rap/jobs', mode: 'copy'

	executor 'local'
	cache true

	cpus = 1
	memory = 4.GB
	time = '10m'

	input:
	tuple val(dataset), val(hash), val(chr), val(jobid)

	output:
	tuple val(dataset), val(hash), val(chr), val(jobid), path("${dataset}-${chr}-${jobid}.json")

	script:
	"""
	dx describe --json ${jobid} > ${dataset}-${chr}-${jobid}.json
	"""

}

/* Clean up failed jobs so that they can be rerun */
process CLEAN {
	tag "${dataset} ${chr}"

	executor 'local'

	cpus = 1
	memory = 4.GB
	time = '10m'

	input:
	tuple val(dataset), val(hash), val(chr), val(jobid), val(json)

	output:
	tuple val(dataset), val(hash), val(chr), val(jobid)

	script:
	"""
	dx unset_properties ${jobid} hash
	"""

}

/* regenie step 2 association testing 
   submit as a run command to RAP
*/
process STEP2 {
	tag "${dataset} ${chr}"

    executor 'local'
	cache true
	
	cpus = 1
	memory = 2.GB
	time = '10m'
	
	input:
	tuple val(dataset), val(hash), val(chr), val(jobtxt), val(traits), val(cluster), val(pheno), val(covar), val(flags), val(pred), val(covar_list), val(covar_cat_list), val(iin)
	
	output:
	tuple val(dataset), val(hash), val(chr), env(JOBID)
	
	script:
	"""
    JOBID=\$(dx run swiss-army-knife \
	--name "regenie ${dataset} ${chr}" \
	--property dataset=${dataset} \
	--property hash=${hash} \
	--property chr=${chr} \
    -iin="/Bulk/Imputation/Imputation from genotype (TOPmed)/ukb21007_c${chr}_b0_v1.bgen" \
    -iin="/Bulk/Imputation/Imputation from genotype (TOPmed)/ukb21007_c${chr}_b0_v1.bgen.bgi" \
    -iin="/Bulk/Imputation/Imputation from genotype (TOPmed)/ukb21007_c${chr}_b0_v1.sample" \
    ${iin.join(' ')} \
    -icmd="curl -O -L https://github.com/rgcgithub/regenie/releases/download/v3.2.6/regenie_v3.2.6.gz_x86_64_Linux_mkl.zip; unzip regenie_v3.2.9.gz_x86_64_Linux_mkl.zip -d /usr/local/bin/; chmod a+x /usr/local/bin/regenie_v3.2.9.gz_x86_64_Linux_mkl; rm regenie*.zip; regenie_v3.2.9.gz_x86_64_Linux_mkl --step 2 ${flags} --bgen ukb21007_c${chr}_b0_v1.bgen --sample ukb21007_c${chr}_b0_v1.sample --phenoFile ${pheno} --covarFile ${covar} --covarColList ${covar_list} --catCovarList ${covar_cat_list} --pred ${pred} --bsize 400 --minMAC 20 --minINFO 0.1 --out step2_${dataset}_chr${chr} --gz --threads 8" \
    --destination="/data/gwas/${dataset}/${hash}" \
    --instance-type "mem1_ssd1_v2_x8" \
    --priority normal \
    --brief --yes)

	"""
}

/* unstage */
process FETCH {
    tag "${dataset}"

    executor 'local'

    cpus = 1
    memory = 4.GB
    time = '30m'

    input:
	tuple val(dataset), val(hash)

    output:
    tuple val(dataset), path("*.regenie.gz"), path("*.log")

    script:
    """
    dx download --recursive "/data/gwas/${dataset}/${hash}/*"
    """

}

/* merge regenie files together for each GWAS.
   sort and rename chromosomes to "chrN"
*/
process MERGE {
	tag "${dataset}-${pheno}"
	
	cpus = 2
	memory = 24.GB
	time = '30m'
	
	input:
	tuple val(dataset), val(pheno), path(regenies)
	val keep
	val remove
	
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
			   col_types=c(EXTRA=col_character())
		)
	) |>
	arrange(CHROM, GENPOS) |>
	mutate(CHROM = if_else(CHROM == 23, true = "X", false = as.character(CHROM)))
	
	write_tsv(regenies, "ukb21007_topmed-${pheno}-${dataset}-${keep}-no${remove}.regenie.gz")
	"""
}

/* collate runtime information */
process RUNTIME {
	publishDir "rap/", mode: 'copy'
	
	cpus = 1
	memory = 4.GB
	time = '10m'
	
	input:
	path(runtime)
	
	output:
	path("rap.tsv")
	
	script:
	"""
	#!Rscript
	library(readr)
	library(dplyr)
	library(lubridate)

	runtime <- read_csv("${runtime}", col_names=c('dataset', 'hash', 'chr', 'jobid', 'start', 'end', 'price'))

	totals <- runtime |>
		mutate(timings=dmilliseconds(end-start)) |>
		group_by(dataset) |>
		summarize(timings=round(sum(timings), 3), price=round(sum(price), 4))

	write_tsv(totals, 'rap.tsv')
	"""
}

/* merge sumstats */
process MANHATTAN {
	tag "${dataset}-${pheno}"
	publishDir "sumstats/${dataset}", mode: 'copy'
	
	cpus = 4
	memory = 48.GB
	time = '2h'
	
	input:
	tuple val(dataset), val(pheno), path(regenie)
	
	output:
	path "*.png"
	
	script:
	"""
	#!Rscript
	library(dplyr)
	library(readr)
	library(fastman)
	
	regenie <- read_tsv("${regenie}") |>
		mutate(P=10^(-LOG10P))
	
	png("${regenie.baseName}.manhattan.png", width=10, height=6, units="in", res=300)
	fastman(regenie, chr = "CHROM", bp = "GENPOS", p = "LOG10P", logp=FALSE, maxP=NULL)
	dev.off()
	
	png("${regenie.baseName}.qq.png", width=6, height=6, units="in", res=300)
	regenie |>
		mutate(p1 = if_else(P == 0, true = .Machine[['double.xmin']], false = P)) |>
		pull(p1) -> p1
	fastqq(p1 = p1, maxP = NULL)
	dev.off()
	"""
}


/* get RSIDs matching the all GWAS CPIDs from the TOPMed sites reference*/
process REGIONS {
	tag "${dataset}"
	
	cpus = 1
	memory = 4.GB
	time = '30m'
	
	input:
	tuple val(dataset), val(pheno), path(regenie), val(dbsnp), path(dbsnp_vcf)
	
	output:
	tuple val(dataset), val(pheno), path(regenie, includeInputs: true), path("${dbsnp_vcf[0].simpleName}-rsids.tsv.gz")
	
	script:
	"""
	# get chrom and position in 0-indexed bed format
	gunzip -c ${regenie} | awk '{OFS="\\t"; if(NR > 1) {print "chr"\$1, \$2-1, \$2}}' > targets.bed
	# block zip and index
	bgzip targets.bed
	tabix targets.bed.gz

	bcftools view --targets-file targets.bed.gz ${dbsnp}.gz |\
	bcftools query -f "%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\n" |\
	split_alts.awk |\
	gzip -c > ${dbsnp_vcf[0].simpleName}-rsids.tsv.gz
	"""
}

/* rename CPID to rsID */
process CPID_RSID {
	tag "${dataset}-${pheno}"
	publishDir "sumstats/${dataset}", mode: 'copy'
	
	cpus = 4
	memory = 48.GB
	time = '2h'
	
	input:
	tuple val(dataset), val(pheno), path(regenie), path(dbsnp)
	
	output:
	path "*.tsv.gz"
	
	script:
	"""
	#!Rscript
	library(dplyr)
	library(readr)
	library(stringr)
	library(fastman)
	
	regenie <- read_tsv("${regenie}", col_types = cols(CHROM=col_character(), GENPOS=col_integer())) 
	dbsnp <- read_tsv("${dbsnp}", comment = "##",
		col_names=c("CHROM", "GENPOS", "RSID", "REF", "ALT"),
		col_types = cols(CHROM = col_character(), GENPOS=col_integer()))
		
	regenie_rsids <- regenie |>
		mutate(CHROM = str_c("chr", CHROM)) |>
		mutate(CPID = if_else(ID == ".", true=str_c(CHROM, ":", GENPOS, ":", ALLELE0, ":", ALLELE1), false=ID)) |>
		left_join(dbsnp, by = c("CHROM", "GENPOS", "ALLELE0"="ALT", "ALLELE1"="REF")) |>
		mutate(ID = coalesce(RSID, CPID)) |>
		select(-RSID, -CPID)

	write_tsv(regenie_rsids, "${regenie.baseName}.tsv.gz")
	"""
}