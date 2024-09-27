nextflow.enable.dsl=2

/*
	Gene-mapping of GWAS meta on the human genome assembly GRCh37 (hg19) from Genome Reference Consortium
*/

// input meta-analysis sumstats
params.meta = "liftover/fixed-*.vcf.gz"

// ld reference PGEN where first phenotype specifies reference population */
params.ref = "reference/ukb_imp_v3.qc_ancestry.{pgen,psam,pvar.zst}"
// gene position reference
params.genes = "reference/glist_ensgid_hg19_v40.txt"
params.names = "reference/glist_ensgid_hg19_v40_symbol_gene_names.txt"

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
	path("pops.txt")

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
	label 'tools'

	cpus = 1
	memory = 16.GB
	time = '30m'

	input: 
	tuple val(pop), val(meta), val(pheno), path(sumstats)
	each neff_pct

	output:
	tuple val(pop), val(meta), val(pheno), path("${sumstats.simpleName}.txt")

	shell:
  	"""
  	echo "SNP\tA1\tA2\tOR\tP\tINFO\tFRQ\tN" > !{sumstats.simpleName}.txt
  	bcftools query \
  	-f "%ID\\t%ALT\\t%REF\\t[%ES]\\t[%LP]\\t[%SI]\\t[%AFCON]\\t[%NE]\\n" \
  	!{sumstats} | awk -v OFS='\t' '{print $1, $2, $3,  exp($4), 10^-($5), $6, $7, $8}' >> !{sumstats.simpleName}.txt
  	"""
}

/* Get genotypes for reference population.
	Rename variants to match input sumstats
*/
process REF_BED {
    tag "${ma.baseName}-${ref}"

    cpus = 1
    memory = 8.GB
    time = '10m'

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

    export PATH=$PATH:/exports/igmm/eddie/GenScotDepression/local/bin/plink2

    plink2 \
    --make-bed \
    --pfile 'vzs' ${ref} \
		--keep-cat-names ${pop} \
		--keep-cat-pheno SuperPop \
		--keep-founders \
    --extract 'bed1' ${ma.baseName}.bed1 \
    --set-all-var-ids @:#:\\\$r:\\\$a \
    --new-id-max-allele-len 500 missing \
    --out ref-cpid \
    --allow-extra-chr \
		--threads ${task.cpus} \
		--memory ${task.memory.bytes.intdiv(1000000)}

    plink2 \
    --make-bed \
    --bfile ref-cpid \
    --update-name ${ma.baseName}.names \
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

    cpus = 8
    memory = 32.GB
    time = '3h'

    input:
    tuple val(pop), val(meta), val(pheno), path(ma), val(ref), path(bedbimbam), path(genelist)

    output:
    tuple val(pop), val(meta), val(pheno), path("${ma.baseName}-${ref}.gene.assoc.mbat"), path("${ma.baseName}-${ref}.gene.snpset.mbat"), path("${ma.baseName}-${ref}.log")

    script:
    """
    /gpfs/igmmfs01/eddie/GenScotDepression/amelia/packages/gcta-1.94.1-linux-kernel-3-x86_64/gcta64 \
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

	publishDir 'maps_hg19', mode: 'copy'

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
	library(readr)
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


