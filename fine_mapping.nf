nextflow.enable.dsl=2

/*
  Fine-mapping of GWAS meta using SuSiEx
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
  // Edit this to read in multiple different ancestries:
  params.ref = "reference/ukb_imp_v3.qc_ancestry.{pgen,psam,pvar.zst}"

} else if (params.build == 'hg38') {
  // files for hg38 build
  // input meta-analysis sumstats
  params.meta = "meta/fixed-*.meta.gz"

  // ld reference PGEN where first phenotype specifies reference population */
  params.ref = "reference/all_hg38.{pgen,psam,pvar.zst}"

}

workflow {

}

