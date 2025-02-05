/* cross cluster genetic correlation estimates */
/* https://github.com/brielin/Popcorn */


params.vcf = "vcf/gwas/GRCh38/*.{csv,json,vcf.gz,vcf.gz.tbi}"
params.reference = "reference/all_hg38.{pgen,psam,pvar.zst}"
params.popcorn = "vendor/Popcorn"

workflow {

    VCF_CH = Channel.fromFilePairs(params.vcf) { it -> it.simpleName }
    REF_CH = Channel.fromFilePairs(params.reference, size: 3) { it -> it.simpleName }

    // reference files
    REF_CLUSTERS_CH = REF_CLUSTERS(REF_CH)
        .splitCsv()
        .map { it -> it[0] }

    REF_QC = QC(REF_CH, REF_CLUSTERS_CH)

    // combine reference files
    REF_REF_CH = REF_QC.combine(REF_QC)
        .filter { it -> it[0] < it[2] }
        .first()
        .view()

    // popscorn scores
    SCORE_CH = COMPUTE(REF_REF_CH)

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
    tuple val(ref), path(pgen)
    each cluster

    output:
    tuple val(cluster), path("${cluster}.*")

    script:
    """
      plink2 \
    --make-bed \
    --pfile 'vzs' ${ref} \
    --keep-cat-names ${cluster} \
    --keep-cat-pheno SuperPop \
    --keep-founders \
    --allow-extra-chr \
    --chr 1-22 \
    --maf 0.05 \
    --rm-dup 'exclude-all' \
    --out ${cluster} \
	--threads ${task.cpus} \
	--memory ${task.memory.bytes.intdiv(1000000)}
    """
}

// compute popcorn scores
process COMPUTE {
    tag "${cluster1}-${cluster2}"

    cpus = 1
    memory = 8.GB

    input:
    tuple val(cluster1), path(bed1), val(cluster2), path(bed2)

    output:
    tuple val(cluster1), val(cluster2), path("scores.txt")

    script:
    """
    popcorn compute -v 1 --bfile1 ${cluster1} --bfile2 ${cluster2} scores.txt
    """
}

