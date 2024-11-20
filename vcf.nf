/**
	Convert sumstats to GWAS VCF format
	
	Inputs:
	  - format/cohort.txt: formatted GWAS sumstats
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
	  - format/cohort.json: sumstats meta data		

**/
import groovy.json.JsonSlurper
import groovy.json.JsonOutput


// Input files: sumstats txt and meta json file
params.sumstats = "format/gwas/*-GRCh38.{txt,json}"

// Output location
params.publish = "vcf"

// reference files
// assembly fasta, fai, and dict 
params.assembly = "reference/Homo_sapiens_assembly38.{fasta,fasta.fai,dict}"
// dbsnp rsID vcf
params.dbsnp = "reference/dbsnp.v153.hg38.vcf.{gz,gz.tbi}"
// chromosome pattenr (if there is one)
params.chr = "chr#"

// gwas2vcf files
params.json = "format/gwas.json"
params.gwas2vcf = "vendor/gwas2vcf"

workflow {

/* 
	Input sumstats
*/

  def jsonSlurper = new JsonSlurper()
	// sumstats and meta data, parse json
	SUMSTATS_CH = Channel
		.fromFilePairs(params.sumstats, size: 2, checkIfExists: true) { it -> it.baseName }
    	.map { it -> it.plus(jsonSlurper.parseText(it[1][0].text)) }
  
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
		.map { it -> jsonSlurper.parseText(it.text) }
	
	GWAS2VCF_CH = Channel
		.fromPath(params.gwas2vcf, type: "dir", checkIfExists: true)

/*
	GWASVCF
*/

	// add genome build to json parameter file
	SUMSTATS_JSON_CH = SUMSTATS_CH
		.combine(JSON_CH)
		.map { it -> it.plus(JsonOutput.toJson(it[3].plus([build: it[2].build]))) }


	// subset sumstats and references file to each chromosome
	// using pattern to get sequence names (usually 'chr#' or '#')
	CHR_CH = Channel.of(1..22, 'X')
		.map { it -> params.chr.replaceFirst('#', it.toString()) }
	
  // qc sumstats
	QC_CH = QC(SUMSTATS_JSON_CH)
  	
	// subset sumstats to chromosome
	QC_CHR_CH = CHR(QC_CH, CHR_CH)
    	.combine(DBSNP_CH)		
		.combine(ASSEMBLY_CH)
		.combine(GWAS2VCF_CH)
    
  
	VCF_CHR_CH = VCF(QC_CHR_CH)
	VCF_TBI_CHR_CH = INDEX(VCF_CHR_CH)
		.map { it -> [it[0], it[2], it[3]] }
		.groupTuple(size: 23)
	
	VCF_CH = CONCAT(VCF_TBI_CHR_CH)

	// additional annotations (AFCAS, AFCON, NE)
	ANNOTATIONS_CH = ANNOTATIONS(QC_CH)

	// harmonised effect size and allele frequency
	MAPPING_CH = MAPPING(VCF_CH)
	
	ANNOTATIONS_MAPPING_CH = ANNOTATIONS_CH
		.join(MAPPING_CH)
		
	ANNOTATIONS_HARMONISED_CH = HARMONISE(ANNOTATIONS_MAPPING_CH)
	
	VCF_ANNOTATIONS_CH = VCF_CH
		.join(ANNOTATIONS_HARMONISED_CH)

	VCF_ANNOTATED_CH = ANNOTATE(VCF_ANNOTATIONS_CH)


}

// QC sumstats for duplicate variants
process QC {
	tag "${dataset}"
  label 'rscript'
	
	cpus = 4
	memory = 32.GB
	time = '1h'
	
	input:
	tuple val(dataset), path(gwas), val(meta), val(json), val(params)
	
	output:
	tuple val(dataset), val(meta), val(params), path("${dataset}-qc.txt")
	
	script:
	"""
	#! /bin/env Rscript
	
	library(readr)
	library(dplyr)
	library(plyranges)
	
	gwas <- read_tsv("${dataset}.txt")
	
	# GRanges object for efficient position lookups
	gwas_gr <- as_granges(gwas, seqnames = chr, start = pos, width = 1)
	
	# count overlaps
	gwas_overlap_counts <- count_overlaps(gwas_gr, gwas_gr)
	
	# find variants with more than 1 overlap where chr+pos+ea/oa are duplicated
	gwas_overlaps <- gwas |> slice(which(gwas_overlap_counts > 1)) |>
		mutate(a1 = pmin(ea, oa), a2 = pmax(ea, oa)) |>
		group_by(chr, pos, a1, a2) |>
		mutate(n = dplyr::n()) |>
		ungroup() |>
		filter(n > 1) |> 
		select(chr, pos, ea, oa)
	
	# remove duplicate positions with an anti-join
	gwas_qc <- gwas |>
		anti_join(gwas_overlaps)
		
	write_tsv(gwas_qc, "${dataset}-qc.txt")
	"""
	
}

// subset to chromosome
process CHR {
	tag "${dataset} ${chr}"
  label 'tools'

	cpus = 1
	memory = 1.GB
	time = '10m'

	input:
	tuple val(dataset), val(meta), val(params), path(gwas)
	each chr

	output:
	tuple val(dataset), val(chr), val(meta), val(params), path("${dataset}-${chr}.txt")

	script:
	"""
	cat ${gwas} | awk 'NR == 1 || \$1 == "${chr}"' > ${dataset}-${chr}.txt
	"""
}

// Convert to VCF
// turn json param string into a file
// run gwas2vcf
process VCF {
	tag "${dataset} ${assembly} ${dbsnp}"
  	label 'gwas2vcf'

	scratch true
	//stageInMode 'copy'
	//stageOutMode 'copy'
  	errorStrategy 'finish'
  

	cpus = 1
	memory = 8.GB
	time = '3h'

	input:
	tuple val(dataset), val(chr), val(meta), val(params), path(gwas), val(dbsnp), path(vcf), val(assembly), path(fasta), path(gwas2vcf)

	output:
	tuple val(dataset), val(meta), path("${dataset}_${chr}.vcf.gz")

	script:
	"""
	cat <<EOF > params.json
	${params}
	EOF
	python ${gwas2vcf}/main.py \
	--data ${gwas} \
	--json params.json \
	--id ${dataset} \
	--ref ${assembly}.fasta \
	--dbsnp ${dbsnp}.gz \
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
	tuple val(dataset), val(dict), path(vcf)

	output:
	tuple val(dataset), val(dict), path(vcf, includeInputs: true), path("${vcf}.tbi")

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
	tag "${dataset}"

	cpus = 1
	memory = 1.GB
	time = '10m'

	input:
	tuple val(dataset), val(dict), val(params), path(sumstats)

	output:
	tuple val(dataset), path("${dataset}-annot.tsv")

	script:
	"""
	cat ${sumstats} | awk 'OFS="\\t" {if(NR == 1) {print "#CHROM", "POS", "REF", "ALT", "BETA", "AFCAS", "AFCON", "NE"} else {print \$1, \$2, \$4, \$3, \$5, \$13, \$14, (4*\$8*\$9)/(\$8+\$9)}}' > ${dataset}-annot.tsv
	"""
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

	headers = c('##FORMAT=<ID=AFCAS,Number=A,Type=Float,Description="Alternate allele frequency in cases">',
			    '##FORMAT=<ID=AFCON,Number=A,Type=Float,Description="Alternate allele frequency in controls">',
	            '##FORMAT=<ID=NE,Number=A,Type=Float,Description="Effective sample size used to estimate genetic effect">')
	write(headers, "header.txt", ncolumns = 1)
	"""
}

// add annotations to GWASVCF
process ANNOTATE {
	tag "${dataset}"
  label 'tools'

	publishDir params.publish, mode: "copy"

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
