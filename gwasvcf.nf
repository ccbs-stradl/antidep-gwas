//nextflow.enable.dsl=2

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
params.meta = "sumstats.csv"
params.fasta = "https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta"
params.fai = "https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai"
params.dict = "https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dict"

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
	
	// run original sumstats through its reformatting script
	FORMAT_CH = FORMAT(SUMSTATS_META_CH)

	// separate sumstats by build
	BUILDS_CH = FORMAT_CH
		.branch {
			GRCh38: it[1].build == "GRCh38"
			GRCh37: it[1].build == "GRCh37"
		}
	
	

}

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
