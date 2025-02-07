/* cross cluster genetic correlation estimates */
/* https://github.com/brielin/Popcorn */


params.vcf = "vcf/gwas/GRCh38/*.{csv,json,vcf.gz,vcf.gz.tbi}"
params.reference = "reference/all_hg38.{pgen,psam,pvar.zst}"
params.exclude = "https://raw.githubusercontent.com/gabraham/flashpca/refs/heads/master/exclusion_regions_hg19.txt"
params.popcorn = "vendor/Popcorn"

import groovy.json.JsonSlurper

workflow {

    def jsonSlurper = new JsonSlurper()
    VCF_CH = Channel.fromFilePairs(params.vcf, size: 4) { it -> it.simpleName }
        .map { it -> it.plus(jsonSlurper.parseText(it[1][1].text)) }
    
    REF_CH = Channel.fromFilePairs(params.reference, size: 3) { it -> it.simpleName }
    EXCLUDE_CH = Channel.fromPath(params.exclude)

    // reference files
    REF_CLUSTERS_CH = REF_CLUSTERS(REF_CH)
        .splitCsv()
        .map { it -> it[0] }

    REF_QC_CH = QC(REF_CH.combine(EXCLUDE_CH), REF_CLUSTERS_CH)

    // combine reference files
    REFS_CH = REF_QC_CH.combine(REF_QC_CH)
        .filter { it -> it[0] <= it[2] }
        .branch { it ->
            eq: it[0] == it[2]
            neq: true
        }
    // remove second set of paths if they are the same
    REF_REF_CH = REFS_CH.eq
        .map { it -> [it[0], it[1], it[2], []] }
        .concat(REFS_CH.neq)
    

    // popscorn scores
    SCORE_CH = COMPUTE(REF_REF_CH)

    format and combine sumstats
    TXT_CH = TXT(VCF_CH)
    TXT_TXT_CH = TXT_CH.combine(TXT_CH)
        .filter { it -> it[0] <= it[3] && it[1] != it[4] }
        .map { it -> [it[0], it[3], it[1], it[2], it[4], it[5]] }

    TXT_SCORE_CH = TXT_TXT_CH
        .combine(SCORE_CH, by: [0, 1])
        .first()

    FIT_CH = FIT(TXT_SCORE_CH)
}

// get reference clusters
process REF_CLUSTERS {
    tag "${ref}"

    cpus = 1

    input:
    tuple val(ref), path(pgen)

    output:
    path("ref.txt")

    script:
    """
    cat ${ref}.psam | awk 'NR > 1 {print \$5}' | sort | uniq > ref.txt
    """
}

// QC reference for each cluster
process QC {
    tag "${cluster}"

    cpus = 1
    memory = 8.GB

    input:
    tuple val(ref), path(pgen), path(exclude)
    each cluster

    output:
    tuple val(cluster), path("${cluster}.*")

    script:
    """
      plink2 \
    --make-bed \
    --set-all-var-ids @:# \
    --pfile 'vzs' ${ref} \
    --keep-cat-names ${cluster} \
    --keep-cat-pheno SuperPop \
    --keep-founders \
    --allow-extra-chr \
    --chr 1-22 \
    --maf 0.05 \
    --rm-dup 'exclude-all' \
    --exclude 'bed1' ${exclude} \
    --out ${cluster} \
	--threads ${task.cpus} \
	--memory ${task.memory.bytes.intdiv(1000000)}
    """
}

// compute popcorn scores
process COMPUTE {
    tag "${cluster1}-${cluster2}"

    cpus = 8
    memory = 16.GB

    input:
    tuple val(cluster1), path(bed1), val(cluster2), path(bed2)

    output:
    tuple val(cluster1), val(cluster2), path("scores.txt")

    script:
    """
    export OMP_NUM_THREADS=${task.cpus}
    popcorn compute -v 1 --bfile1 ${cluster1} --bfile2 ${cluster2} scores.txt
    """
}

// format sumstats
process TXT {
  tag "${dataset}"
  label 'tools'
  
  cpus = 1
  memory = 1.GB
  
  input:
  tuple val(dataset), path(vcf), val(info)
  
  output:
  tuple val(info.cluster), val(dataset), path("*.txt")
  
  shell:
  '''
  echo "rsid\ta1\ta2\tN\tbeta\tSE\taf" > !{dataset}.txt
  bcftools query \
  -f "%CHROM %POS %ALT %REF [%NE] [%ES] [%SE] [%AFCON]\\n" \
  !{dataset}.vcf.gz | awk -v OFS='\t' '{sub(/chr/, "", $1); print $1":"$2, $3, $4, $5, $6, $7, 1-$8}' >> !{dataset}.txt
  '''
}

// fit popcorn corrleations
process FIT {
    tag "${dataset1}_${dataset2}"

    publishDir "models/popcorn", mode: "copy"

    cpus = 8
    memory = 16.GB

    input:
    tuple val(cluster1), val(cluster2), val(dataset1), path(sumstats1), val(dataset2), path(sumstats2), path(scores)

    output:
    tuple val(cluster1), val(cluster2), val(dataset1), val(dataset2), path("${dataset1}--${dataset2}.correlations.txt")

    script:
    """
    export OMP_NUM_THREADS=${task.cpus}
    popcorn fit -v 1 --cfile ${scores} --sfile1 ${sumstats1} --sfile2 ${sumstats2} ${dataset1}--${dataset2}.correlations.txt
    """
}