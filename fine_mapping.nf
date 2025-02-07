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
nextflow log condescending_swartz -f hash,process,tag

*/
import groovy.json.JsonSlurper
import groovy.json.JsonOutput
import java.nio.file.Paths

// MAKE_BFILE INPUTS:
  params.pfile = '/exports/igmm/eddie/GenScotDepression/data/ukb/genetics/impv3_pgen/ukb_imp_v3.qc' // prefix of the pfiles to pass to --pfile in plink 
  params.ancestry_ids = 'reference/ukb-ld-ref_ancestry.id' // Output from script make-pgen.sh
  params.chr = 3 // testing pipeline on chromosome I know there are results for, so i dont have to deal with NULL clumping file

  // effective sample size QC parameter
  params.neff_pct = 0.8 // gwas is filtered by effectve sample size, threshold is 0.8

workflow {

/*
  MAKE_BFILE process
  make bim, bed, fam files for each ancestry/cluster
  this only needs to happen for chromosomes that will be fine mapped, 
  is there a way to get this from the sumstats?
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

  // MA INPUTS:
    // load in gwas meta sumstats that have been lifted over to hg19 and match ancestries in CLUSTER_CH

    META_CH = CLUSTER_CH 
      .map { cluster -> "vcf/meta/GRCh37/antidep-2501-fixed-N06A-${cluster}.vcf.gz" } // edit this line so path is not hard coded
      .map { pathStr -> Paths.get(pathStr) } // can this be replaced with a nextflow command?
      .map { it -> [it.simpleName.split("-"), it] }
      .map { it -> [it[0][4], it[0][2], it[0][3], it[1]] } // pull out sections of file name: [4] = ancestry/cluster, [2] = fixed, [3] = N06A.

	// QC parameters
	NEFF_CH = Channel.of(params.neff_pct)


  // format sumstats to .ma
	MA_CH = MA(META_CH, NEFF_CH)

/*
  Extract effective sample size for each meta-analysis sumstats
*/
  // Effective sample size from the meta analysis for each ancestry/cluster
  // Extract neff values for each ancestry


  // NEFF INPUTS:
   // json contains metadata about the meta analysis gwas, eg. pheno, ancestry, sample sizes, effective sample size

  def jsonSlurper = new JsonSlurper()

  NEFF_TOTAL_CH = CLUSTER_CH
    .map { cluster -> "vcf/meta/GRCh37/antidep-2501-fixed-N06A-${cluster}.json" } // edit this line so path is not hard coded
    .map { path -> 
        def jsonFile = new File(path)
        jsonSlurper.parseText(jsonFile.text)
    }
    .map { json_obj -> 
        [json_obj.cluster, json_obj.neff] 
    }


/*
  Combine ancestry/cluster, with MA (path to sumstats), and total NEFF
  eg. channel will have the first element as a list of ancestries,
  second element will be path to sumstats (in same order as ancestries in first element)
  third element will be effective sample size (in same order as ancestries in first element)
  fourth element will be bfile prefix
  these will all be string values that are comma separated (with no spaces around commas, else susiex throws an error)

  #########
  Example code: 
  JOINED_CH = Channel.of(
    ['EUR', 'PATH1', '100'],
    ['AFR', 'PATH2', '200'],
    ['SAS', 'PATH3', '300'] 
  )

  JOINED_CH
    .collate(3)
    .transpose()
    .map { it.join(',') } 
    .collect()
    .view()

  Outputs: ['EUR,AFR,SAS', 'PATH1,PATH2,PATH3', '100,200,300']
  #########

*/

  JOINED_CH = MA_CH
    .join(NEFF_TOTAL_CH) // joined based on ancestry/cluster value (first element of channel)
    .map { it -> [it[0], it[3], it[4]] } // keep only ancestry/cluster value, path to MA and neff total.
    .join(BFILE_CH.map { it -> [it[0], "${file(it[2][0]).parent}/${file(it[2][0]).baseName}" ] }) // keep only ancestry/cluster value and bfile (without file extension)
    .collate(4)
    .transpose()
    .map { it.join(',') } 
    .collect()

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

  // Combine all the ancestries clumping results into one channel,
  // as we need all ancestries in the R script together to 
  // find overlapping regions between ancestries.


  CLUMP_CLUSTER_CH = CLUMP_CH
    .collect(flat : false)

  CLUMP_POST_CH = CLUMP_POST(CLUMP_CLUSTER_CH)

  JOINED_SUSIEX_CH = CLUMP_POST_CH
    .map { it ->  [it[0], file(it[1]).text.trim()] } // read value of chr.txt
    .combine(JOINED_CH)

/*
  SUSIEX process
  run SuSiEx on each region  
*/

  SUSIEX_CH = SUSIEX(JOINED_SUSIEX_CH)

/*
  SUSIEX_POST process
  explore results using susiexR package
*/

  SUSIEX_PROCESSED_CH = SUSIEX_CH
                          .map { it -> [ file(it[0][0]).parent , it[4], it[5], it[6] ]}

  SUSIEX_PROCESSED_CH.view()

  SUSIEX_POST_CH = SUSIEX_POST(SUSIEX_PROCESSED_CH)

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

    Rscript ${baseDir}/fine_mapping_ma_process.R "${sumstats.simpleName}.ma"
	  """

}

process CLUMP {
  tag "${chr}.${ma}"
  
  cpus = 1
  memory = 8.GB
  time = '10m'

  input:
    tuple val(cluster), val(meta), val(pheno), path(ma), val(bfile_prefix), path(bfiles), val(chr)

  output:
    tuple val(cluster), val(meta), val(pheno), path(ma), path("${ma.baseName}.${chr}.clumps"), path("${ma.baseName}.${chr}.log"), val(chr)
 
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


process CLUMP_POST {
  // would be good to tag with chromosome, but nervous to change the input to this process as it was fidly to get right

  label 'analysis'

  cpus = 1
  memory = 32.GB
  time = '30m'

  input:
    val(nested_clumps)

  output:
    tuple path("*.finemapRegions"), path("chr.txt")

  script:
  """
  Rscript ${baseDir}/fine_mapping_clump_post.R "${nested_clumps}"
  """

} // Warning: if edits are made to this script nextflow won't know to re-run it if output is already cached


process SUSIEX {
  tag "chr:${chr}"

  cpus = 4
  memory = 32.GB
  time = '30m'

  publishDir "fineMapping/output", mode: "copy"

  input:
    tuple path(finemapRegions), val(chr), val(ancestries), val(maPaths), val(neff), val(bfile)

  output:
    tuple path("*.log"), path("*.cs"), path("*.snp"), path("*.summary"), val(maPaths), val(ancestries), val(chr)

  script:
  """

  CLUMP_RANGES_FILE=${finemapRegions}

  tail -n +2 "\${CLUMP_RANGES_FILE}" | while IFS=\$'\\t' read -r line; do
    # Get CHR, BP_START, BP_END by using awk to match column names indices
    # Really important we put these cols in the correct order in the R scripts making the "loci" object
    CHR=\$(echo "\$line" | awk '{print \$1+0}')
    BP_START=\$(echo "\$line" | awk '{print \$2+0}')
    BP_END=\$(echo "\$line" | awk '{print \$3+0}')

    echo "Processing CHR: \$CHR, BP_START: \$BP_START, BP_END: \$BP_END"

    SuSiEx \
     --sst_file=${maPaths} \
     --n_gwas=${neff} \
     --ref_file=${bfile} \
     --ld_file=${ancestries} \
     --out_dir=. \
     --out_name=SuSiEx.${ancestries}.output.cs95_\${CHR}:\${BP_START}:\${BP_END} \
     --level=0.95 \
     --pval_thresh=1e-5 \
     --chr=\$CHR \
     --bp=\$BP_START,\$BP_END \
     --maf=0.005 \
     --snp_col=1,1,1 \
     --chr_col=9,9,9 \
     --bp_col=10,10,10 \
     --a1_col=2,2,2 \
     --a2_col=3,3,3 \
     --eff_col=5,5,5 \
     --se_col=6,6,6 \
     --pval_col=7,7,7 \
     --mult-step=True \
     --plink=plink \
     --keep-ambig=True |& tee SuSiEx.${ancestries}.output.cs95_\${CHR}:\${BP_START}:\${BP_END}.log

  done 

  """

} 


process SUSIEX_POST {
  tag "chr:${chr}"
  label 'rscript'

  cpus = 2
  memory = 32.GB
  time = '30m'

  publishDir "fineMapping/plots", mode: "copy"

  input:
    tuple val(susiexPath), val(sumstatsPath), val(ancestries), val(chr)

  output:
    path("*.png"), optional: true

  script:
  """
  #!Rscript

  # Ensure all R libraries are installed prior to nextflow execution for version of R loaded
  #if (!require("susiexr", quietly = TRUE)) { devtools::install_github("ameliaes/susiexr") }
  library(cowplot)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(susiexR)

  # ---- Format SuSiEx results:
  path_to_susiex_results <- "${susiexPath}"

  ancestries <- str_split("${ancestries}", ",")[[1]]

  results <- format_results(path_to_susiex_results, ancestries = ancestries)

  # Check that each processed file type has the same number of fine mapped regions
  nrow(results\$summary) == length(results\$cs) && length(results\$cs) == length(results\$snp)

  # ---- Plot the relationship between:
  # CS_LENGTH - number of SNPs in the credible set
  # CS_PURITY - purity of the credible set
  # MAX_PIP - Maximum posterior inclusion probability (PIP) in the credible set.

  png("length_purity_maxPIP.png", width = 1000, height = 600, res = 150)
    print(plotPurityPIP(results\$summary))
  dev.off()

  # ---- Plot the probability the top SNP in the credible set is causal in each ancestry
  # ----------------
  png("POST-HOC_PROB_POP.png", width = 1000, height = 800, res = 150)
    print(plotAncestryCausal(results\$summary, ancestries = ancestries))
  dev.off()



  """

}


/*
  



*/
