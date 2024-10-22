nextflow.enable.dsl=2

/*
	Gene-mapping of GWAS meta
*/

// genotyype build parameters
// first define which options are available
def valid_builds = ['hg19', 'hg38']

// default to hg19 if not provided
params.build = params.build ?: 'hg19'

// validate that the provided build is one of the valid options
if (!(params.build in valid_builds)) {
    error "Invalid value for --build. Please specify one of: ${valid_builds.join(', ')}"
}

// assign labels for process MA based on build:
def processMALabel = (params.build == 'hg19') ? 'tools' : (params.build == 'hg38') ? 'rscript' : null

def mapsDuplicatesDir = (params.build == 'hg19') ? 'maps_hg19_duplicates' : (params.build == 'hg38') ? 'maps_hg38_duplicates' : null

def mapsDir = (params.build == 'hg19') ? 'maps_hg19' : (params.build == 'hg38') ? 'maps_hg38' : null

// input files dependent on genotype build type:
if (params.build == 'hg19') {
	// files for hg19 build
	// input meta-analysis sumstats
	params.meta = "liftover/fixed-*.vcf.gz"

	// ld reference PGEN where first phenotype specifies reference population */
	params.ref = "reference/ukb_imp_v3.qc_ancestry.{pgen,psam,pvar.zst}"

	// gene position reference
	params.genes = "reference/glist_ensgid_hg19_v40.txt"
	params.names = "reference/glist_ensgid_hg19_v40_symbol_gene_names.txt"

} else if (params.build == 'hg38') {
	// files for hg38 build
	// input meta-analysis sumstats
	params.meta = "meta/fixed-*.meta.gz"

	// ld reference PGEN where first phenotype specifies reference population */
	params.ref = "reference/all_hg38.{pgen,psam,pvar.zst}"

	// gene position reference
	params.genes = "reference/glist_ensgid_hg38_v40.txt"
	params.names = "reference/glist_ensgid_hg38_v40_symbol_gene_names.txt"
}

// effective sample size QC parameter
params.neff_pct = 0.8

workflow {

	// sumstats from fixed effects meta 
	META_CH = Channel.fromPath(params.meta)
		.map { it -> [it.simpleName.split("-"), it] }
		.map { it -> [it[0][2], it[0][0], it[0][1], it[1]] }

	// reference population PBFILEs 
	REF_CH = Channel.fromFilePairs(params.ref, size: 3)

	// populations represented in the reference file
	POP_CH = POPS(REF_CH)
		.splitCsv(header: ['pop'])
		.map { it -> it.pop }

	// gene name and position reference
	GENES_CH = Channel.fromPath(params.genes)
	GENENAMES_CH = Channel.fromPath(params.names)

	// QC parameters
	NEFF_CH = Channel.of(params.neff_pct)

	// format sumstats to .ma
	MA_CH = MA(META_CH, NEFF_CH)

	// retain sumstats to analyze that are in the reference sample
	META_POP_CH = POP_CH
		.cross(MA_CH)
		.map { it -> it[1] }
		.combine(REF_CH)

	META_POP_BED_CH = REF_BED(META_POP_CH)
		.combine(GENES_CH)

	MBAT_CH = MBAT(META_POP_BED_CH)
		.combine(GENENAMES_CH)

	MBAT_RESULTS_CH = MBATTED(MBAT_CH)
}

/* Get reference populations from first phenotype in reference genotypes files */
process POPS {
	tag "${ref}"

	cpus = 1
	memory = 4.GB
	time = '10m'

	executor 'local'

	input:
	tuple val(ref), path(pgensamvar)

	output:
	path("pops.txt"), emit: pops_txt

	script:
	"""
	cat ${ref}.psam | awk '{if(NR > 1) {print \$5}}' | sort | uniq > pops.txt
	"""

}

/* Reformat to .ma for GCTA input 
	exclude SNPs with low relative NEFF
*/
process MA {
	tag "${sumstats}"
	label processMALabel

	cpus = 1
	memory = 16.GB
	time = '30m'

	input: 
	tuple val(pop), val(meta), val(pheno), path(sumstats)
	each neff_pct

	output:
	tuple val(pop), val(meta), val(pheno), path("${sumstats.simpleName}.ma")

	script:
	if (params.build == 'hg19') {
		shell:
	  	"""
		echo -e "SNP\tA1\tA2\tfreq\tBETA\tSE\tP\tN\tCHR\tBP\tNE" > ${sumstats.simpleName}.ma
		bcftools query \
	    -f "%ID\\t%ALT\\t%REF\\t[%AFCON]\\t[%ES]\\t[%SE]\\t[%LP]\\t[%SS]\\t%CHROM\\t%POS\\t[%NE]" \
	    ${sumstats} | awk -v OFS='\\t' -v neff_threshold=!{params.neff_pct} '\$11 >= neff_threshold * \$11 {print \$1, \$2, \$3, \$4, \$5, \$6, 10^-(\$7), \$8, \$9, \$10}' >> ${sumstats.simpleName}.ma
	  	"""
	} else if (params.build == 'hg38') {
		"""
		#!Rscript
		if( !require('readr') ){
		install.packages('readr', repos = "https://cloud.r-project.org/")
		}
		library(readr)

		if( !require('dplyr') ){
		install.packages('dplyr', repos = "https://cloud.r-project.org/")
		}
		library(dplyr)

		if( !require('stringr') ){
		install.packages('stringr', repos = "https://cloud.r-project.org/")
		}
		library(stringr)

		sumstats <- read_tsv("${sumstats}")

		ma <- sumstats |>
			filter(NEFF >= ${neff_pct} * NEFF) |>
			transmute(SNP, A1, A2, freq = AFCON, BETA = log(OR), SE, P, N = NEFF,
					str_remove(CHR, "chr"), BP)
		
		write_tsv(ma, "${sumstats.simpleName}.ma")
		"""
	}
}

/* Get genotypes for reference population.
	Rename variants to match input sumstats
*/
process REF_BED {
	label "plink2"
    tag "${ma.baseName}-${ref}"

    cpus = 1
    memory = 8.GB
    time = '10m'

	publishDir mapsDuplicatesDir, mode: 'copy', pattern: '*.duplicates'

    input:
    tuple val(pop), val(meta), val(pheno), path(ma), val(ref), path(pgen)

    output:
    tuple val(pop), val(meta), val(pheno), path(ma, includeInputs: true), val("ref-${ref}"), path("ref-${ref}.{bed,bim,fam}")

    script:
    """
    # chr/pos to extract
    cat ${ma} | awk 'NR > 1 {print \$9, \$10, \$10}' > ${ma.baseName}.bed1

    # rename CPIDs to rsID
    cat ${ma} | awk 'NR > 1 {print \$9":"\$10":"\$3":"\$2, \$1}' > ${ma.baseName}.names

	# Save duplicate CPIDs
    awk '{print \$1}' ${ma.baseName}.names | sort |uniq -d > ${ma.baseName}.duplicates

	# Remove duplicate CPIDs from *.names, if there are duplicates
	if [ ! -s ${ma.baseName}.duplicates ]; then
    	# If empty, just copy the original .names file to .noDups.names
    	cp ${ma.baseName}.names ${ma.baseName}.noDups.names
	else
    	# Remove duplicate CPIDs from *.names
    	awk 'NR==FNR {duplicates[\$1]; next} !(\$1 in duplicates)' ${ma.baseName}.duplicates ${ma.baseName}.names > ${ma.baseName}.noDups.names
	fi

	plink2 \
    --make-bed \
    --pfile 'vzs' ${ref} \
		--keep-cat-names ${pop} \
		--keep-cat-pheno SuperPop \
		--keep-founders \
    --extract 'bed1' ${ma.baseName}.bed1 \
    --set-all-var-ids @:#:\\\$r:\\\$a \
    --new-id-max-allele-len 700 missing \
    --out ref-cpid \
    --allow-extra-chr \
		--threads ${task.cpus} \
		--memory ${task.memory.bytes.intdiv(1000000)}

    plink2 \
    --make-bed \
    --bfile ref-cpid \
    --update-name ${ma.baseName}.noDups.names \
    --out ref-${ref} \
		--threads ${task.cpus} \
		--memory ${task.memory.bytes.intdiv(1000000)}
    """
}

/* Gene-based testing: mBAT-combo
   https://yanglab.westlake.edu.cn/software/gcta/#mBAT-combo
*/
process MBAT {
    tag "${ma.baseName}-${ref}"
    label "gcta64"

    cpus = 8
    memory = 32.GB
    time = '3h'

    input:
    tuple val(pop), val(meta), val(pheno), path(ma), val(ref), path(bedbimbam), path(genelist)

    output:
    tuple val(pop), val(meta), val(pheno), path("${ma.baseName}-${ref}.gene.assoc.mbat"), path("${ma.baseName}-${ref}.gene.snpset.mbat"), path("${ma.baseName}-${ref}.log")

    script:
    """
     gcta64_v1.94 \
		--bfile ${ref} \
		--mBAT-combo ${ma} \
		--mBAT-gene-list ${genelist} \
		--mBAT-print-all-p \
		--mBAT-svd-gamma 0.9 \
		--mBAT-wind 50 \
		--diff-freq 0.2 \
		--fastBAT-ld-cutoff 0.9 \
		--mBAT-write-snpset \
		--out ${ma.baseName}-${ref} \
		--threads ${task.cpus}
    """
}

/* Get gene names, FDR-correct 
*/
process MBATTED {
	tag "${mbat.simpleName}"
	label "rscript"

	publishDir mapsDir , mode: 'copy'

	cpus = 1
	memory = 4.GB
	time = '10m'

	input:
	tuple val(pop), val(meta), val(pheno), path(mbat), path(snpset), path(log), path(genes)

	output:
	tuple val(pop), val(meta), val(pheno), path("${mbat.simpleName}.mbat.tsv")

	script:
	"""
	#!Rscript
	if( !require('readr') ){
	install.packages('readr', repos = "https://cloud.r-project.org/")
	}
	library(readr)

	if( !require('dplyr') ){
	install.packages('dplyr', repos = "https://cloud.r-project.org/")
	}
	library(dplyr)

	mbat <- read_tsv("${mbat}")
	genes <- read_tsv("${genes}")

	mbat_genes <- mbat |>
		left_join(select(genes, ensgid, gene_name), by=c("Gene" = "ensgid")) |>
		filter(P_mBATcombo <= 2.7e-6) |>
		select(Gene, gene_name,  everything())

	write_tsv(mbat_genes, "${mbat.simpleName}.mbat.tsv")
	"""
}


