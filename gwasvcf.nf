nextflow.preview.dsl=2

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

params.sumstats = "sumstats/*.{gz,sh}"
params.fasta = "https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta"
params.fai = "https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai"
params.dict = "https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dict"

workflow {

	// sumstats and reformating script
	SUMSTATS_CH = Channel
	.fromFilePairs(params.sumstats, size: 2, checkIfExists: true)
	
	FORMAT_CH = FORMAT(SUMSTATS_CH)
	.view()

}

process FORMAT {
	tag "${cohort}"

	executor = 'local'

	cpus = 1
	memory = 1.GB
	time = '10m'

	input:
	tuple val(cohort), path(sumstats)

	output:
	tuple val(cohort), path("${cohort}.txt")

	shell:
	"""
	sh ${cohort}.sh ${cohort}.gz ${cohort}.txt
	"""
}
