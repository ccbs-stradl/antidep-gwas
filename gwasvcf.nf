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
params.dbsnp = "reference/dbsnp.v153.hg38.vcf.{gz,gz.tbi}"
// hg19 to hg38 liftover chain
params.chain = "reference/hg19ToHg38.over.chain.gz"

// gwas2vcf files
params.json = "sumstats/gwas.json"
params.gwas2vcf = "vendor/gwas2vcf"

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
		.map {it -> ["${it[0][1].cohort}-${it[0][1].pheno}-${it[0][1].cluster}", it[0][1], it[0][0], it[1][1]] }

/*
	Reference files
*/

	ASSEMBLY_CH = Channel
		.fromFilePairs(params.assembly, size: 3, checkIfExists: true)

	DBSNP_CH = Channel
		.fromFilePairs(params.dbsnp, size: 2, checkIfExists: true)

	CHAIN_CH = Channel
		.fromPath(params.chain, checkIfExists: true)

/*
	gwas2vcf files
*/
	JSON_CH = Channel
		.fromPath(params.json, checkIfExists: true)
	
	GWAS2VCF_CH = Channel
		.fromPath(params.gwas2vcf, type: "dir", checkIfExists: true)

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
			[it[0], [cohort: it[1].cohort, pheno: it[1].pheno, filekey: it[1].filekey, build: "GRCh38", cluster: it[1].cluster], it[2]]
		}

/*
	GWASVCF
*/

	// merge b38 sumstats with lifted-over b37 sumstats
	// tack on reference dbsnp vcf
	// and gwas2vcf files
	DATA_CH = BUILDS_CH.GRCh38
		.concat(BUILDS_37_TO_38_CH)
		.combine(DBSNP_CH)		
	// subset sumstats and references file to each chromosome
	CHR_CH = Channel.of(1..22, 'X')
		.map { it -> "chr${it}" }
	
	DATA_CHR_CH = CHR(DATA_CH, CHR_CH)
		.combine(ASSEMBLY_CH)
		.combine(JSON_CH)
		.combine(GWAS2VCF_CH)

	VCF_CHR_CH = VCF(DATA_CHR_CH)
	VCF_TBI_CHR_CH = INDEX(VCF_CHR_CH)
		.map { it -> [it[0], it[2], it[3]] }
		.groupTuple()
		
	VCF_CH = CONCAT(VCF_TBI_CHR_CH)

	// additional annotations (AFCAS, AFCON, NE)
	ANNOTATIONS_CH = ANNOTATIONS(DATA_CH)
	
	VCF_ANNOTATIONS_CH = VCF_CH
		.join(ANNOTATIONS_CH)

	VCF_ANNOTATED_CH = ANNOTATE(VCF_ANNOTATIONS_CH)

}

// Reformat sumstats for GWASVCF input
process FORMAT {
	tag "${dataset}"

	cpus = 1
	memory = 1.GB
	time = '10m'

	input:
	tuple val(dataset), val(dict), val(filekey), path(sumstats)

	output:
	tuple val(dataset), val(dict), path("${dataset}.txt")

	shell:
	"""
	sh ${filekey}.sh ${filekey}.gz ${dataset}.txt
	"""
}

// get rsIDs for GWAS variants from dbSNP
process SELECT {
	tag "${dataset} ${dbsnp}"
	
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'

	cpus = 1
	memory =32.GB
	time = '4h'

	input:
	tuple val(dataset), val(dict), path(gwas), val(dbsnp), path(dbsnp_vcf)

	output:
	tuple val(dataset), val(dict), path(gwas, includeInputs: true), path("${dataset}-rsids.tsv")

	script:
	"""
	# get list of SNPs from the GWAS
	cat ${gwas} | awk '{if(NR > 1) {print \$10}}' > ${dataset}.list

	# extract variants from the dbSNP reference based on rsID
	gatk SelectVariants \
	-V ${dbsnp}.gz \
	--keep-ids ${dataset}.list \
	-O ${dataset}-select.vcf.gz

	# get chrom, pos, and id for selected rsIDs
	bcftools query \
	-f "%CHROM\\t%POS\\t%ID\n" \
	${dataset}-select.vcf.gz > ${dataset}-rsids.tsv
	"""
}


// split sumstats into parts that can/cannot be lifted by rsID
// output a cpid list to be lifted in bed format
process PRELIFT {
	tag "${dataset}"

	cpus = 1
	memory =16.GB
	time = '4h'

	input:
	tuple val(dataset), val(dict), path(gwas), path(rsid)

	output:
	tuple val(dataset), val(dict), path(gwas, includeInputs: true), path(rsid, includeInputs: true), path("${dataset}-cpid.bed")

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

	write_tsv(bed_lift_by_cpid, "${dataset}-cpid.bed", col_names=FALSE)
	"""
}

// liftover cpids using chain file
process LIFTOVER {
	tag "${dataset} ${chain}"

	cpus = 1
	memory =4.GB
	time = '1h'

	input:
	tuple val(dataset), val(dict), path(gwas), path(rsid), path(cpid), path(chain)

	output:
	tuple val(dataset), val(dict), path(gwas), path(rsid), path("${dataset}-cpid-lifted.bed"), path("${dataset}-cpid-unlifted.bed")

	script:
	"""
	liftOver \
	${cpid} \
	${chain} \
	${dataset}-cpid-lifted.bed \
	${dataset}-cpid-unlifted.bed
	"""
}


// update sumstats with liftover files
process LIFT {
	tag "${dataset}"

	cpus = 1
	memory =16.GB
	time = '1h'

	input:
	tuple val(dataset), val(dict), path(gwas), path(rsid), path(lifted), path(unlifted)

	output:
	tuple val(dataset), val(dict), path("${gwas.simpleName}-lifted.txt")

	script:
	"""
	#!/bin/env Rscript

	library(readr)
	library(dplyr)
	library(stringr)

	gwas <- read_tsv("${gwas}")
	rsid_liftover <- read_tsv("${rsid}", col_names = c("chr", "pos", "snp"))
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

// subset to chromosome
process CHR {
	tag "${dataset} ${dbsnp} ${chr}"

	cpus = 1
	memory =4.GB
	time = '3h'

	input:
	tuple val(dataset), val(dict), path(gwas), val(dbsnp), path(vcf)
	each chr

	output:
	tuple val(dataset), val(dict), val(chr), path("${dataset}-${chr}.txt"), val("${vcf[0].simpleName}_${dataset}_${chr}"), path("${vcf[0].simpleName}_${dataset}_${chr}.vcf.{gz,gz.tbi}")

	script:
	"""
	cat ${gwas} | awk 'NF == 1 || \$1 == "${chr}"' > ${dataset}-${chr}.txt
	cat ${dataset}-${chr}.txt | awk '{OFS="\\t"; {print \$1, \$2-1, \$2}}' > ${dataset}_${chr}.bed
	bgzip ${dataset}_${chr}.bed
	tabix ${dataset}_${chr}.bed.gz
	
	bcftools view --targets-file ${dataset}_${chr}.bed.gz -o ${vcf[0].simpleName}_${dataset}_${chr}.vcf.gz ${dbsnp}.gz
	tabix ${vcf[0].simpleName}_${dataset}_${chr}.vcf.gz
	"""
}

// Convert to VCF
process VCF {
	tag "${dataset} ${assembly} ${dbsnp}"

	container = "docker://mrcieu/gwas2vcf:latest"
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'

	cpus = 1
	memory =16.GB
	time = '48h'

	input:
	tuple val(dataset), val(dict), val(chr), path(gwas), val(dbsnp), path(vcf), val(assembly), path(fasta), path(json), path(gwas2vcf)

	output:
	tuple val(dataset), val(dict), path("${dataset}_${chr}.vcf.gz")

	script:
	"""
	python ${gwas2vcf}/main.py \
	--data ${gwas} \
	--json ${json} \
	--id ${dataset} \
	--ref ${assembly}.fasta \
	--dbsnp ${dbsnp}.vcf.gz \
	--out ./${dataset}_${chr}.vcf.gz
	"""
}

// Index VCF
process INDEX {
	tag "${vcf.simpleName}"

	cpus = 1
	memory =4.GB
	time = '30m'

	input:
	tuple val(dataset), val(dict), path(vcf)

	output:
	tuple val(dataset), val(dict), path(vcf, includeInputs: true), path("${vcf}.tbi")

	script:
	"""
	tabix ${vcf}
	"""
}

process CONCAT {
	tag "${dataset}"

	cpus = 1
	memory = 4.GB
	time = '1h'

	input:
	tuple val(dataset), path(vcf), path(tbi)

	output:
	tuple val(dataset), path("${dataset}-base.vcf.gz"), path("${dataset}-base.vcf.gz.tbi")

	script:
	"""
	bcftools concat --output ${dataset}-base.vcf.gz ${vcf}
	tabix ${dataset}-base.vcf.gz
	"""
}

process ANNOTATIONS {
	tag "${dataset}"

	cpus = 1
	memory = 1.GB
	time = '10m'

	input:
	tuple val(dataset), val(dict), val(sumstats), val(dbsnp), path(vcf)

	output:
	tuple val(dataset), val("${dataset}-annot"), path("${dataset}-annot.{gz,gz.tbi}")

	script:
	"""
	cat ${sumstats} | awk 'OFS="\\t" {if(NR == 1) {print "#CHROM", "POS", "AFCAS", "AFCON", "NE"} else {print \$1, \$2, \$13, \$14, (4*\$8*\$9)/(\$8+\$9)}}' | bgzip -c > ${dataset}-annot.gz
	tabix -s1 -b2 -e2 ${dataset}-annot.gz
	"""
}

process ANNOTATE {
	tag "${dataset}"

	publishDir "vcf", mode: "copy"

	cpus = 1
	memory = 4.GB
	time = '1h'

	input:
	tuple val(dataset), path(vcf), path(tbi), val(annot), path(annotations)

	output:
	tuple val(dataset), path("${dataset}.vcf.gz"), path("${dataset}.vcf.gz.tbi")

	script:
	"""
	echo -e '##FORMAT=<ID=AFCAS,Number=1,Type=Integer,Description="Alternate allele frequency in cases">' > hdr.txt
	echo -e '##FORMAT=<ID=AFCON,Number=1,Type=Integer,Description="Alternate allele frequency in controls">' >> hdr.txt
	echo -e '##FORMAT=<ID=NE,Number=1,Type=Integer,Description="Effective sample size used to estimate genetic effect">' >> hdr.txt
	
	bcftools annotate -s ${dataset} -a ${annot}.gz -h hdr.txt -c CHROM,POS,FORMAT/AFCAS,FORMAT/AFCON,FORMAT/NE -o ${dataset}.vcf.gz ${vcf}

	tabix ${dataset}.vcf.gz
	"""
}