nextflow.enable.dsl=2

/*
  Fine-mapping of GWAS meta using SuSiEx (hg19 only)

  Summary of processes and inputs for each process:
    * make bim, bed, fam files for each ancestry:
      - pfile 
      - IDs where for each ancestry (output from make-pgen.sh)

    * process meta sumstats into correct format:
      - meta-analysed sumstats: meta/metaset-method-phenotype-cluster.csv (output from format_meta.nf)

    * calculate effective sample size for each meta-analysis sumstats:
      - format/meta/GRCh38/antidep-2408-fixed-N06A-{ANCESTRY}.json (output from format_meta.nf)

    * clumping to identify regions for fine mapping:
      - processed sumstats (from previous process in current script)
      - bfile, per chr and ancestry, (output from make-pgen_hg19.sh)

    * load clumping results into R and determine region boundaries:
      - output from clumping process, per chr and ancestry (from previous process in current script)

    * run SuSiEx on each region:
      - output from R determining region boundaries (from previous process in current script)
      - processed sumstats (from process above)

    * explore results:
      - use susiexR package in R to format and plot output from SuSiEx


// nextflow log gives the dir for each process, 
nextflow log fabulous_ampere -f hash,process 
nextflow log fabulous_ampere -f hash,process,tag

*/
import groovy.json.JsonSlurper

// MAKE_BFILE INPUTS:
  params.pfile = '/exports/igmm/eddie/GenScotDepression/data/ukb/genetics/impv3_pgen/ukb_imp_v3.qc' // prefix of the pfiles to pass to --pfile in plink 
  params.ancestry_ids = 'reference/ukb-ld-ref_ancestry.id' // Output from script make-pgen.sh
  params.chr = 21 // testing pipeline on a small chromosome

// MA INPUTS:
  params.meta = "liftover/fixed-*.vcf.gz" // load in gwas meta sumstats that have been lifted over to hg19

  // effective sample size QC parameter
  params.neff_pct = 0.8 // gwas is filtered by effectve sample size, threshold is 0.8

// NEFF INPUTS:
  params.neff_total = 'format/meta/GRCh38/antidep-2408-fixed-N06A-EUR.json' // this json contains metadata about the meta analysis gwas, eg. pheno, ancestry, sample sizes, effective sample size

workflow {

/*
  MAKE_BFILE process
  make bim, bed, fam files for each ancestry/cluster
*/
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

/*
  MA process
  process meta sumstats into correct format
*/
  // Get sumstats from fixed effects meta 
  // Make a channel with contents: val(cluster), val(meta), val(pheno), path(sumstats)
    META_CH = Channel.fromPath(params.meta)
      .map { it -> [it.simpleName.split("-"), it] }
      .map { it -> [it[0][2], it[0][0], it[0][1], it[1]] }

	// QC parameters
	NEFF_CH = Channel.of(params.neff_pct)

  // format sumstats to .ma
	MA_CH = MA(META_CH, NEFF_CH)

/*
  NEFF_TOTAL process
  calculate effective sample size for each meta-analysis sumstats
*/
  // Effective sample size from the meta analysis for each ancestry/cluster
  // Extract neff values for each ancestry

  def jsonSlurper = new JsonSlurper()

  NEFF_TOTAL_CH = Channel.fromPath(params.neff_total)
    .map { file -> jsonSlurper.parseText(file.text) }
  
  NEFF_TOTAL(NEFF_TOTAL_CH)

/*
  CLUMP process
  clumping to identify regions for fine mapping
*/
    // Creates channel with:
    // val(cluster), val(meta), val(pheno), path("${sumstats.simpleName}.ma"), path(pgen), path(psam), path(pvar)
    MA_BFILE_CH = MA_CH
      .join(BFILE_CH)

  CLUMP_CH = CLUMP(MA_BFILE_CH)

/*
  CLUMP_POST process
  load clumping results into R and determine region boundaries
*/

  CLUMP_CH.view()

  // CLUMP_POST(CLUMP_CH)

/*
  SUSIEX process
  run SuSiEx on each region  
*/



/*
  SUSIEX_POST process
  explore results using susiexR package
*/


}

process MAKE_BFILE {
  tag "$cluster:chr$chr"

  cpus = 1
  memory = 8.GB
  time = '1h'

  input:
    tuple val(pfile), path(ancestry_ids), val(cluster), val(chr)

  output:
    tuple val(cluster), val("ukb_imp_v3.qc.geno02.mind02_${cluster}_${chr}"), path("ukb_imp_v3.qc.geno02.mind02_${cluster}_${chr}.*"), val(chr)
 
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
	tuple val(cluster), val(meta), val(pheno), path(sumstats)
	each neff_pct

	output:
	tuple val(cluster), val(meta), val(pheno), path("${sumstats.simpleName}.ma")

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
    tuple val(cluster), val(meta), val(pheno), path(ma), val(bfile_prefix), path(bfiles), val(chr)

  output:
    tuple val(cluster), val(meta), val(pheno), path(ma), path("${ma.baseName}.${chr}.clumps"), path("${ma.baseName}.${chr}.log")
 
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
    --bfile ${bfile_prefix} \
    --out ${ma.baseName}.${chr} \
    --threads ${task.cpus} \
	  --memory ${task.memory.bytes.intdiv(1000000)}
  
  # Check if .clumps file exists; if not, create the file with content "NULL"
  # This is needed otherwise nextflow exits with an error because plink has not made a .clumps file (because there are no clumps in the gwas for that chromosome)
    if [ ! -f "${ma.baseName}.${chr}.clumps" ]; then
        echo "NULL" > "${ma.baseName}.${chr}.clumps"
    fi

  """
}

process CHECK_NEFF {
  tag "json"
  
  cpus = 1
  memory = 8.GB
  time = '10m'

  input:
  val(jsonData)

  output:
  stdout

  script:
  """
  echo "The effective sample size is ${jsonData.neff}"

  # Check output from this by checking the directory:
  # nextflow log fabulous_ampere -f hash,process,tag
  # cat ../../ad/work/2c/1bc695013dc56de81a4fd7fa5fd453/.command.out 
  """

}


process CLUMP_POST {
  tag "${result}"
  label 'analysis'

  cpus = 1
  memory = 32.GB
  time = '30m'

  input:
    tuple val(cluster), val(meta), val(pheno), path(ma), path(clumps), path(log)

  output:


  script:
  """
  #!Rscript
  library(data.table)
  library(dplyr)
  library(plyranges) # reduce_ranges

  # Read in clumps results
  clumped_data <- fread("${clumps}")

  # TODO:
  # If clumped_data == NULL then skip file processing
  #
  #
  #
  #
  #

  loci <- clumped_data %>%
    dplyr::select(CHR= `#CHROM`, POS, SNP=ID)

  hg19 <- genome_info('hg19')
  # pull out lengths manually since seqnames uses "chrN" instead of "N"
  hg19_chr_lengths <- as_tibble(hg19) |> slice(1:23) |> pull(width)

  # add genome info for autosomes and X
  grng <- loci %>%
            arrange(CHR) %>%
            as_granges(seqnames = CHR,
                          start = POS,
                          end = POS) 

  seqlevels(grng) <- as.character(1:23)
  seqlengths(grng) <- hg19_chr_lengths
  grng <- set_genome_info(grng, 'hg19', is_circular = rep(FALSE, 23))

  # stretch each region, by 00 kb upstream and downstream, then trim back to position boundaries
  grng_streched <- stretch(anchor_center(grng), 200000) %>% trim()

  # TO DO
  # Example code to get overlapping regions for multiple ancestries
  #
  #
  #
  #
  #
  #
  #

  # Reduce ranges to collapse overlapping or nearby regions
  grng_reduced <- grng_streched %>%
    reduce_ranges(min.gapwidth = 5000)  %>% # What should this be set to?
    as_tibble() %>%
    mutate(WIDTH = end - start + 1) %>%
    dplyr::select(CHR = seqnames, BP_START = start, BP_END = end, WIDTH)


  # Save the table of min and max BP positions for SuSiEx
  write.table(grng_reduced, "test/fixed-N06A-EUR.human_g1k_v37.neff08_noZero.clumpRanges",
  row.names = F, quote = F, sep = "\t")
  """

}


/*

process SUSIEX {

# output:
SuSiEx.${ANCESTRY_NAMES_CONCAT}.output.cs95_${CHR}:${BP_START}:${BP_END}.log

# publishDir = fineMapping

# ensure SuSiEx and plink have been copied to bin folder

    SuSiEx \ 
      --sst_file= ${SUMSTAT_FILES_ORDERED_BY_ANCESTRY} \
      --n_gwas= ${NEFF_TOTAL_ORDERED_BY_ANCESTRY} \
      --ref_file= ${BFILE_PREFIX_ORDERED_BY_ANCESTRY, PER CHROMOSOME} \
      --ld_file= ${PATH_TO_TMP_DIR_ORDERED_BY_ANCESTRY, PER CHROMOSOME} \
      --out_dir=./fineMapping/results \
      --out_name=SuSiEx.${ANCESTRY_NAMES_CONCAT}.output.cs95_${CHR}:${BP_START}:${BP_END} \
      --level=0.95 \
      --pval_thresh=1e-5 \
      --chr=$CHR \
      --bp= $BP_START,$BP_END \
      --maf=0.005 \
      --snp_col=1,1,1 \ REP EACH COL NUMBER BY LENGTH OF ANCESTRIES
      --chr_col=9,9,9 \
      --bp_col=10,10,10 \
      --a1_col=2,2,2 \
      --a2_col=3,3,3 \
      --eff_col=5,5,5 \
      --se_col=6,6,6 \
      --pval_col=7,7,7 \
      --mult-step=True \
      --plink=plink \
      --keep-ambig=True |& tee SuSiEx.${ANCESTRY_NAMES_CONCAT}.output.cs95_${CHR}:${BP_START}:${BP_END}.log


}

process SUSIEX_POST {

}

*/