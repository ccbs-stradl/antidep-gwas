nextflow.enable.dsl=2

/*
  Fine-mapping of GWAS meta using SuSiEx (hg19 only)

  Summary of processes and inputs for each process:
    * make bim, bed, fam files for each ancestry:
      - pfile 
      - IDs where for each ancestry (output from make-pgen.sh)
    * calculate effective sample size for each meta-analysis sumstats:
      - format/meta/GRCh38/antidep-2408-fixed-N06A-{ANCESTRY}.json (output from format_meta.nf)
    * process meta sumstats into correct format:
      - meta-analysed sumstats: meta/metaset-method-phenotype-cluster.csv (output from format_meta.nf)
    * clumping to identify regions for fine mapping:
      - processed sumstats (from previous process in current script)
      - pfile, per chr and ancestry, (output from make-pgen_hg19.sh)
    * load clumping results into R and determine region boundaries:
      - output from clumping process, per chr and ancestry (from previous process in current script)
    * run SuSiEx on each region:
      - output from R determining region boundaries (from previous process in current script)
      - processed sumstats (from process above)
    * explore results:
      - use susiexR package in R to format and plot output from SuSiEx

*/

// MAKE_BFILE INPUTS:
params.pfile = '/exports/igmm/eddie/GenScotDepression/data/ukb/genetics/impv3_pgen/ukb_imp_v3.qc'
params.ancestry_ids = 'reference/ukb-ld-ref_ancestry.id'
params.chr = 21

// MA INPUTS:
  params.meta = "liftover/fixed-*.vcf.gz"

  // effective sample size QC parameter
  params.neff_pct = 0.8

// NEFF INPUTS:
  params.neff_total = 'format/meta/GRCh38/antidep-2408-fixed-N06A-EUR.json'

// CLUMP INPUTS:
  params.clump_file_prefix = 'reference/ukb_imp_v3.qc.geno02_'

workflow {

  // Make reference files in bfile format
  PFILE_CH = Channel.of(params.pfile)          
  IDS_CH = Channel.fromPath(params.ancestry_ids)
  CHR_CH = Channel.of(params.chr)    

    // Create CLUSTER_CH to extract unique ancestries from the 3rd column of IDS_CH
    CLUSTER_CH = IDS_CH
      .splitCsv(header: false, sep: ' ') 
      .map { it[2] }                     
      .unique()    

  REF_CH = PFILE_CH
    .combine(IDS_CH)
    .combine(CLUSTER_CH)
    .combine(CHR_CH)

  BFILE_CH = MAKE_BFILE(REF_CH)

  // Get sumstats from fixed effects meta 
  // Make a channel with contents: val(pop), val(meta), val(pheno), path(sumstats)
    META_CH = Channel.fromPath(params.meta)
      .map { it -> [it.simpleName.split("-"), it] }
      .map { it -> [it[0][2], it[0][0], it[0][1], it[1]] }

	// QC parameters
	NEFF_CH = Channel.of(params.neff_pct)

  // format sumstats to .ma
	MA_CH = MA(META_CH, NEFF_CH)

  // Effective sample size from the meta analysis for each ancestry
  // Extract neff values for each ancestry
  // NEFF_TOTAL_CH = 
  Channel.fromPath(params.neff_total)
    .splitJson()
    .filter { it -> it.key == "neff" }
    .map { it.value }
    .view()

  // Clumping to identify regions for fine mapping
      // - processed sumstats (from previous process in current script)
      // - bfile (output from MAKE_BFILE)

    // Creates channel with:
    // val(pop), val(meta), val(pheno), path("${sumstats.simpleName}.ma"), path(pgen), path(psam), path(pvar)
    MA_BFILE_CH = MA_CH
      .join(BFILE_CH)
      .map { pop, meta, pheno, ma, bfile_paths, chr ->
        def prefix = "${bfile_paths[0].getParent()}/${bfile_paths[0].getBaseName()}" 
        return tuple(pop, meta, pheno, ma, prefix, chr)
      }

  // MA_BFILE_CH.view()

  CLUMP_CH = CLUMP(MA_BFILE_CH)

}

process MAKE_BFILE {
  tag "$cluster:chr$chr"

  cpus = 1
  memory = 8.GB
  time = '1h'

  input:
    tuple val(pfile), path(ancestry_ids), val(cluster), val(chr)

  output:
    tuple val(cluster), path("ukb_imp_v3.qc.geno02.mind02_${cluster}_${chr}.*"), val(chr)
 
  script:
  """
  plink2 \
    --pfile  ${pfile} \
    --keep-col-match ${ancestry_ids} ${cluster} \
    --chr ${chr} \
    --geno 0.02 \
    --maf 0.005 \
    --make-pgen 'vzs' \
    --out ukb_imp_v3.qc.geno02_${cluster}_${chr} \
    --threads ${task.cpus} \
	  --memory ${task.memory.bytes.intdiv(1000000)}

    plink2 \
    --pfile ukb_imp_v3.qc.geno02_${cluster}_${chr} 'vzs' \
    --mind 0.02 \
    --make-bed \
    --out ukb_imp_v3.qc.geno02.mind02_${cluster}_${chr} \
    --threads ${task.cpus} \
	  --memory ${task.memory.bytes.intdiv(1000000)}
  """

}

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
	tuple val(pop), val(meta), val(pheno), path("${sumstats.simpleName}.ma")

	script:
	  """
		echo -e "SNP\tA1\tA2\tfreq\tBETA\tSE\tP\tN\tCHR\tBP\tNE" > ${sumstats.simpleName}.ma
		bcftools query \
	    -f "%ID\\t%ALT\\t%REF\\t[%AFCON]\\t[%ES]\\t[%SE]\\t[%LP]\\t[%SS]\\t%CHROM\\t%POS\\t[%NE]" \
	    ${sumstats} | awk -v OFS='\\t' -v neff_threshold=!{params.neff_pct} '\$11 >= neff_threshold * \$11 {print \$1, \$2, \$3, \$4, \$5, \$6, 10^-(\$7), \$8, \$9, \$10, \$11}' >> ${sumstats.simpleName}.ma
  """

}


process CLUMP {
  tag "${ma}"
  
  cpus = 1
  memory = 8.GB
  time = '10m'

  input:
    tuple val(pop), val(meta), val(pheno), path(ma), val(bfile), val(chr)

  output:
    tuple val(pop), val(meta), val(pheno), path(ma), path("${ma.baseName}.${chr}.clumps"), path("${ma.baseName}.${chr}.log")
 
  script:
  """
  plink2 \
    --clump ${ma} \
    --clump-id-field SNP \
    --clump-p-field P \
    --clump-p1 5e-8 \
    --clump-p2 0.05 \
    --clump-r2 0.1 \
    --clump-kb 1000 \
    --bfile ${bfile} \
    --out ${ma.baseName}.${chr} \
    --threads ${task.cpus} \
	  --memory ${task.memory.bytes.intdiv(1000000)}
  
  # Check if .clumps file exists; if not, create the file with content "NULL"
    if [ ! -f "${ma.baseName}.${chr}.clumps" ]; then
        echo "NULL" > "${ma.baseName}.${chr}.clumps"
    fi

  """
}

