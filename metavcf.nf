/* 
	Convert meta-analysis to GWAS VCF format
	
	Inputs:
		- meta/fixed-*.meta.gz: fixed-effects meta-analysis sumstats
		
*/

// Input files: sumstats gz, shell scripts, and meta data csv
params.sumstats = "meta/*.gz"

// reference files
// assembly fasta, fai, and dict 
params.assembly = "reference/Homo_sapiens_assembly38.{fasta,fasta.fai,dict}"
// dbsnp rsID vcf
params.dbsnp = "reference/dbsnp.v153.hg38.vcf.{gz,gz.tbi}"

// gwas2vcf files
params.json = "sumstats/gwas.json"
params.gwas2vcf = "vendor/gwas2vcf"

workflow {

/* 
	Input sumstats
*/

	// meta-analysed sumstats
	SUMSTATS_CH = Channel
		.fromPath(params.sumstats)


/*
	Reference files
*/

	ASSEMBLY_CH = Channel
		.fromFilePairs(params.assembly, size: 3, checkIfExists: true)

	DBSNP_CH = Channel
		.fromFilePairs(params.dbsnp, size: 2, checkIfExists: true)

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
	
	// run original sumstats through a reformatting
	FORMAT_CH = FORMAT(SUMSTATS_CH)	
/*
	GWASVCF
*/

	// subset sumstats and references file to each chromosome
	CHR_CH = Channel.of(1..22, 'X')
		.map { it -> "chr${it}" }

	DATA_CH = FORMAT_CH
		.combine(DBSNP_CH)
		
	DATA_CHR_CH = CHR(DATA_CH, CHR_CH)
		.combine(ASSEMBLY_CH)
		.combine(JSON_CH)
		.combine(GWAS2VCF_CH)

	VCF_CHR_CH = VCF(DATA_CHR_CH)
	VCF_TBI_CHR_CH = INDEX(VCF_CHR_CH)
		.groupTuple(size: 23)
		
	VCF_CH = CONCAT(VCF_TBI_CHR_CH)

	// additional annotations (AFCAS, AFCON, NE)
	ANNOTATIONS_CH = ANNOTATIONS(DATA_CH)

	// harmonised effect size and allele frequency
	MAPPING_CH = MAPPING(VCF_CH)
	
	ANNOTATIONS_MAPPING_CH = ANNOTATIONS_CH
		.join(MAPPING_CH)
	     
	ANNOTATIONS_HARMONISED_CH = HARMONISE(ANNOTATIONS_MAPPING_CH)
	
	VCF_ANNOTATIONS_CH = VCF_CH
		.join(ANNOTATIONS_HARMONISED_CH)

	VCF_ANNOTATED_CH = ANNOTATE(VCF_ANNOTATIONS_CH)
	

}

/* Reformat sumstats for GWASVCF input */
process FORMAT {
	tag "${sumstats.simpleName}"

	cpus = 1
	memory = 1.GB
	time = '10m'

	input:
	path(sumstats)

	output:
	tuple val(sumstats.simpleName), path("*.txt")

	shell:
	'''
	
	# sumstats columns
	#  1  CHR
    #  2  BP
    #  3  SNP
    #  4  A1
    #  5  A2
    #  6  studies
    #  7  OR
    #  8  SE
    #  9  P
    # 10  OR_R
    # 11  SE_R
    # 12  P_R
    # 13  Q
    # 14  I
    # 15  INFO
    # 16  AFCAS
    # 17  AFCON
    # 18  NCAS
    # 19  NCON
    # 20  NEFF
    # 21  NTOT


	gunzip -c !{sumstats} | awk 'OFS = "\t" {if(NR == 1) {print "chr", "pos", "ea", "oa", "beta", "se", "pval", "ncase", "ncontrol", "snp", "eaf", "imp_info", "eaf_case", "eaf_control", "neff"} else {print $1, $2, $4, $5, log($7), $8, $9, $18, $19, $3, $17, $15, $16, $17, $20}}' > !{sumstats.simpleName}.txt

	# output columns
	# chr
	# pos
	# ea
	# oa
	# beta
	# se
	# pval
	# ncase
	# ncontrol
	# snp
	# eaf
	# imp_info
	# eaf_case
	# eaf_control
	# neff
	'''
}


// subset to chromosome
process CHR {
	tag "${dataset} ${dbsnp} ${chr}"
  	label 'tools'

	cpus = 1
	memory =4.GB
	time = '3h'

	input:
	tuple val(dataset), path(gwas), val(dbsnp), path(vcf)
	each chr

	output:
	tuple val(dataset), val(chr), path("${dataset}_${chr}.txt"), val("${vcf[0].simpleName}_${dataset}_${chr}"), path("${vcf[0].simpleName}_${dataset}_${chr}.vcf.{gz,gz.tbi}")

	script:
	"""
	cat ${gwas} | awk 'NR == 1 || \$1 == "${chr}"' > ${dataset}_${chr}.txt
	cat ${dataset}_${chr}.txt | tail -n +2 | awk '{OFS="\\t"; {print \$1, \$2-1, \$2}}' > ${dataset}_${chr}.bed
	bgzip ${dataset}_${chr}.bed
	tabix ${dataset}_${chr}.bed.gz
	
	bcftools view --targets-file ${dataset}_${chr}.bed.gz -o ${vcf[0].simpleName}_${dataset}_${chr}.vcf.gz ${dbsnp}.gz
	tabix ${vcf[0].simpleName}_${dataset}_${chr}.vcf.gz
	"""
}

// Convert to VCF
process VCF {
	tag "${dataset} ${assembly} ${dbsnp}"
  	label 'gwas2vcf'

	//scratch true
	//stageInMode 'copy'
	//stageOutMode 'copy'
  errorStrategy 'finish'
  

	cpus = 1
	memory = 8.GB
	time = '1h'

	input:
	tuple val(dataset), val(chr), path(gwas), val(dbsnp), path(vcf), val(assembly), path(fasta), path(json), path(gwas2vcf)

	output:
	tuple val(dataset), path("${dataset}_${chr}.vcf.gz")

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
  label 'tools'

	cpus = 1
	memory =4.GB
	time = '30m'

	input:
	tuple val(dataset), path(vcf)

	output:
	tuple val(dataset), path(vcf, includeInputs: true), path("${vcf}.tbi")

	script:
	"""
	tabix ${vcf}
	"""
}

// Concatenate per-chromosome VCFs together
process CONCAT {
	tag "${dataset}"
  label 'tools'

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

// make annotations file
process ANNOTATIONS {
	tag "${sumstats.baseName}"

	cpus = 1
	memory = 1.GB
	time = '10m'

	input:
	tuple val(dataset), path(sumstats), val(dbsnp), path(vcf)

	output:
	tuple val(dataset), path("${dataset}-annot.tsv")

	shell:
	'''	
	cat !{sumstats} | awk 'OFS="\\t" {if(NR == 1) {print "#CHROM", "POS", "REF", "ALT", "BETA", "AFCAS", "AFCON", "NE"} else {print $1, $2, $4, $3, $5, $13, $14, $15}}' > !{dataset}-annot.tsv
	'''
}

// get harmonised effect size and allele frequency
process MAPPING {
	tag "${dataset}"
  label 'tools'

	cpus = 1
	memory = 1.GB
	time = '10m'

	input:
	tuple val(dataset), path(vcf), path(tbi)

	output:
	tuple val(dataset), path("${dataset}-mapping.tsv")

	script:
	"""
	echo -e "#CHROM\\tPOS\\tREF\\tALT\\tES\\tAF" > ${dataset}-mapping.tsv
	bcftools query -f "%CHROM\\t%POS\\t%REF\\t%ALT\\t[%ES]\\t[%AF]\\n" ${vcf} >> ${dataset}-mapping.tsv
	"""
}

// annotations harmonised with what is in the VCF file
// output tsv that can be tabix'd and a list of annotation column names
process HARMONISE {
	tag "${dataset}"
  label 'rscript'
	
	cpus = 4
	memory = 32.GB
	time = '30m'
	
	input:
	tuple val(dataset), path(annotations), path(mapping)
	
	output:
	tuple val(dataset), path("${dataset}-annot-mapped.tsv"), path("${dataset}-annot-mapped.cols"), path("header.txt")
	
	script:
	"""
	#!/bin/env Rscript
	
	library(readr)
	library(dplyr)
	library(stringr)
	
	# effect size and allele frequency from the harmonised VCF file
	mapping <- read_tsv("${mapping}")
	# annotations to add
	annot <- read_tsv("${annotations}")
	
	mapping_annot_ref_alt <- mapping |>
		inner_join(annot, by = c("#CHROM" = "#CHROM", "POS"  = "POS", "REF" = "REF", "ALT" = "ALT"))
		
	mapping_annot_alt_ref <- mapping |>
		inner_join(annot, by = c("#CHROM" = "#CHROM", "POS"  = "POS", "ALT" = "REF", "REF" = "ALT"))

	# compare ES to BETA to check if the effect allele has been swapped
	mapped_annot <- bind_rows(mapping_annot_ref_alt, mapping_annot_alt_ref) |>
		mutate(flip = sign(ES) != sign(BETA)) |>
		transmute(CHROM = `#CHROM`, POS, REF, ALT,
		          `FORMAT/AFCAS` = if_else(flip, true = 1 - AFCAS, false = AFCAS),
				      `FORMAT/AFCON` = if_else(flip, true = 1 - AFCON, false = AFCON),
				      `FORMAT/NE` = NE) |>
		arrange(CHROM, POS)

	write_tsv(mapped_annot, "${dataset}-annot-mapped.tsv", col_names = FALSE)
	write(names(mapped_annot), "${dataset}-annot-mapped.cols", ncolumns = 1)

	headers = c('##FORMAT=<ID=AFCAS,Number=1,Type=Float,Description="Alternate allele frequency in cases">',
			        '##FORMAT=<ID=AFCON,Number=1,Type=Float,Description="Alternate allele frequency in controls">',
	            '##FORMAT=<ID=NE,Number=1,Type=Float,Description="Effective sample size used to estimate genetic effect">')
	write(headers, "header.txt", ncolumns = 1)
	"""
}

// add annotations to GWASVCF
process ANNOTATE {
	tag "${dataset}"
  label 'tools'

	publishDir "metavcf", mode: "copy"

	cpus = 1
	memory = 4.GB
	time = '1h'

	input:
	tuple val(dataset), path(vcf), path(tbi), path(annotations), path(columns), path(headers)

	output:
	tuple val(dataset), path("${dataset}.vcf.gz"), path("${dataset}.vcf.gz.tbi")

	script:
	"""
	# compress and index the annotations
	bgzip -c ${annotations} > ${annotations}.gz
	tabix -s 1 -b 2 -e 2 ${annotations}.gz

	# add annotations to the VCF
	bcftools annotate \
	--samples ${dataset} \
	--annotations ${annotations}.gz \
	--columns-file ${columns} \
	--header-lines ${headers} \
	--output ${dataset}.vcf.gz \
	${vcf}

	tabix ${dataset}.vcf.gz
	"""
}
