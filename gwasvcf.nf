nextflow.enable.dsl=2

/* 
	Convert sumstats to GWAS VCF format
	
	Inputs:
		- sumstats/cohort.gz: original GWAS sumstats
		- sumstats/cohort.sh: shell script to reformat sumstats to the following columns:		
			- chr
			- pos
			- ea
			- oa
			- beta
			- se
			- pval
			- ncase
			- ncontrol
			- snp
			- eaf
			- imp_info
			- eaf_case
			- eaf_control
		
*/

// Input files: sumstats gz, shell scripts, and meta data csv
params.sumstats = "sumstats/*.{gz,sh}"
params.meta = "sumstats.csv"

// reference files
// assembly fasta, fai, and dict 
params.assembly = "reference/Homo_sapiens_assembly38.{fasta,fasta.fai,dict}"
// dbsnp rsID vcf
params.dbsnp = "reference/All_20180418.vcf.{gz,gz.tbi}"
// hg19 to hg38 liftover chain
params.chain = "reference/hg19ToHg38.over.chain.gz"

workflow {

/* 
	Input sumstats
*/

	// sumstats and reformating script, keyed to basename
	SUMSTATS_CH = Channel
		.fromFilePairs(params.sumstats, size: 2, checkIfExists: true)

	// sumstats meta data, keyed to filekey column
	META_CH = Channel
		.fromPath(params.meta)
		.splitCsv(header: true)
		.map { it -> [it.filekey, it] }

	// merge sumstats and meta data information
	SUMSTATS_META_CH = META_CH
		.cross(SUMSTATS_CH)
		.map {it -> [it[0][1].cohort, it[0][1], it[0][0], it[1][1]] }

/*
	Reference files
*/

	ASSEMBLY_CH = Channel
		.fromFilePairs(params.assembly, size: 3)

	DBSNP_CH = Channel
		.fromFilePairs(params.dbsnp, size: 2) { it -> it.simpleName }

	CHAIN_CH = Channel
		.fromPath(params.chain)

/*
	Process sumstats
*/
	
	// run original sumstats through its reformatting script
	FORMAT_CH = FORMAT(SUMSTATS_META_CH)

	// separate sumstats by build
	BUILDS_CH = FORMAT_CH
		.branch {
			GRCh38: it[1].build == "GRCh38"
			GRCh37: it[1].build == "GRCh37"
		}


/*
	Liftover for GRCh37 
*/
	// merge GRCh37 with GRCh38 dbSNP reference
	BUILDS_37_CH = BUILDS_CH.GRCh37
		.combine(DBSNP_CH)

	// get GRCh38 CHR/POS for GWAS SNPs based on rsID
	BUILDS_37_SELECT_CH = SELECT(BUILDS_37_CH)

	// splits sumstats into that can/can't be lifted by rsID
	BUILDS_37_PRELIFT_CH = PRELIFT(BUILDS_37_SELECT_CH)
		.combine(CHAIN_CH)

	// perform liftover of cpid list
	BUILDS_37_LIFTED_CH = LIFTOVER(BUILDS_37_PRELIFT_CH)
	
	// update sumstats to new build
	BUILDS_37_TO_38_CH = LIFT(BUILDS_37_LIFTED_CH)
		.map { it -> 
			[it[0], [cohort: it[1].cohort, filekey: it[1].filekey, build: "GRCh38", cluster: it[1].cluster], it[2]]
		}

	DATA_CH = BUILDS_CH.GRCh38
		.concat(BUILDS_37_TO_38_CH)
		.combine(ASSEMBLY_CH)
		.combine(DBSNP_CH)

/*
	GWASVCF
*/

	//VCF_CH = VCF(FASTA_CH, FAI_CH, DICT_CH)

}

// Reformat sumstats for GWASVCF input
process FORMAT {
	tag "${cohort}"

	executor = 'local'

	cpus = 1
	memory = 1.GB
	time = '10m'

	input:
	tuple val(cohort), val(dict), val(filekey), path(sumstats)

	output:
	tuple val(cohort), val(dict), path("${cohort}.txt")

	shell:
	"""
	sh ${filekey}.sh ${filekey}.gz ${cohort}.txt
	"""
}

// get rsIDs for GWAS variants from dbSNP
process SELECT {
	tag "${cohort} ${dbsnp}"

	cpus = 1
	memory =4.GB
	time = '4h'

	input:
	tuple val(cohort), val(dict), path(gwas), val(dbsnp), path(dbsnp_vcf)

	output:
	tuple val(cohort), val(dict), path(gwas, includeInputs: true), path("${cohort}-rsids.tsv")

	script:
	"""
	# get list of SNPs from the GWAS
	cat ${gwas} | awk '{if(NR > 1) {print \$10}}' > ${cohort}.list

	# extract variants from the dbSNP reference based on rsID
	gatk SelectVariants \
	-V ${dbsnp}.vcf.gz \
	--keep-ids ${cohort}.list \
	-O ${cohort}-select.vcf.gz

	# get chrom, pos, and id for selected rsIDs
	bcftools query \
	-f "%CHROM\\t%POS\\t%ID\n" \
	${cohort}-select.vcf.gz > ${cohort}-rsids.tsv
	"""
}


// split sumstats into parts that can/cannot be lifted by rsID
// output a cpid list to be lifted in bed format
process PRELIFT {
	tag "${cohort}"

	cpus = 1
	memory =16.GB
	time = '4h'

	input:
	tuple val(cohort), val(dict), path(gwas), path(rsid)

	output:
	tuple val(cohort), val(dict), path(gwas, includeInputs: true), path(rsid, includeInputs: true), path("${cohort}-cpid.bed")

	script:
	"""
	#! /bin/env Rscript

	library(readr)
	library(dplyr)

	gwas <- read_tsv("${gwas}")
	rsids <- read_tsv("${rsid}", col_names=c('chr', 'pos', 'snp'))

	gwas_lift_by_rsid <- gwas |> filter(snp %in% pull(rsids, snp))
	gwas_lift_by_cpid <- gwas |> filter(!snp %in% pull(rsids, snp))

	bed_lift_by_cpid <- gwas_lift_by_cpid |>
		transmute(chr, pos-1, pos, snp)

	write_tsv(bed_lift_by_cpid, "${cohort}-cpid.bed", col_names=FALSE)
	"""
}

// liftover cpids using chain file
process LIFTOVER {
	tag "${cohort} ${chain}"

	cpus = 1
	memory =4.GB
	time = '1h'

	input:
	tuple val(cohort), val(dict), path(gwas), path(rsid), path(cpid), path(chain)

	output:
	tuple val(cohort), val(dict), path(gwas), path(rsid), path("${cohort}-cpid-lifted.bed"), path("${cohort}-cpid-unlifted.bed")

	script:
	"""
	liftOver \
	${cpid} \
	${chain} \
	${cohort}-cpid-lifted.bed \
	${cohort}-cpid-unlifted.bed
	"""
}


// update sumstats with liftover files
process LIFT {
	tag "${cohort}"

	cpus = 1
	memory =4.GB
	time = '1h'

	input:
	tuple val(cohort), val(dict), path(gwas), path(rsid), path(lifted), path(unlifted)

	output:
	tuple val(cohort), val(dict), path("${gwas.simpleName}-lifted.txt")

	script:
	"""
	#! /bin/env Rscript

	library(readr)
	library(dplyr)
	library(stringr)

	gwas <- read_tsv("${gwas}")
	rsid_liftover <- read_tsv("${rsid}", col_names = c("chr", "pos", "snp")) |>
		mutate(chr=paste0("chr", chr))
	cpid_liftover <- read_tsv("${lifted}", col_names = c("chr", "start", "end", "snp"))

	# merge by snp and replace chr/pos with liftover positions
	gwas_lift_with_rsid <- gwas |>
		inner_join(rsid_liftover, by = "snp", suffix = c(".b37", ".b38"))
		
	gwas_lifted_with_rsid <- gwas_lift_with_rsid |>
		select(chr = chr.b38, pos = pos.b38, everything(), -chr.b37, -pos.b37)

	# merge by snp, replace chr/pos with liftover, positions, keep rsid names and
	# update cpid names (chr.b38:pos.b38:oa:ea)
	gwas_lift_with_cpid <- gwas |>
		inner_join(cpid_liftover, by = "snp", suffix = c(".b37", ".b38"))
		
	gwas_lifted_with_cpid <- gwas_lift_with_cpid |>
		select(chr = chr.b38, pos = end, everything(), -chr.b37, -pos, -start) |>
		mutate(snp = if_else(str_detect(snp, "rs"), 
							 true = snp,
							 false = paste(chr, pos, oa, ea, sep=":")))

	gwas_lifted <- bind_rows(gwas_lifted_with_rsid, gwas_lifted_with_cpid) |>
		mutate(chrom=if_else(chr == "chrX",
					 true = 23,
					 false = as.numeric(str_replace(chr, "chr", "")))) |>
		filter(!is.na(chrom)) |>
		arrange(chrom, pos) |>
		select(-chrom)

	write_tsv(gwas_lifted, "${gwas.simpleName}-lifted.txt")
	"""
}

process VCF {

	cpus = 1
	memory =16.GB
	time = '4h'

	input:
	path(fasta)
	path(fai)
	path(dict)

	output:
	path("*.vcf")

	script:
	"""
	ls ${fasta}
	ls ${fai}
	ls ${dict}
	"""
}