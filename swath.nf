nextflow.enable.dsl=2

/*
    GWASVCF files in Hail
*/

params.vcf = "vcf/*.vcf.gz"
params.name = "swath"

workflow {
    VCF_CH = Channel
        .fromPath(params.vcf)

    NAME_CH = Channel
        .of(params.name)

    // Index the VCFs
    VCF_TBI_CH = INDEX(VCF_CH)
        .multiMap { it ->
            vcf: it[0]
            tbi: it[1]    
        }

    // Merge all VCFs together
    VCF_ALL_CH = VCF_TBI_CH
        .vcf
        .collect()

     TBI_ALL_CH = VCF_TBI_CH
        .tbi
        .collect()

    MERGE_CH = MERGE(VCF_ALL_CH, TBI_ALL_CH, NAME_CH)

    // Import VCF as a Hail MatrixTable
    MT_CH = MT(MERGE_CH)
}

process INDEX {
    tag "${vcf.simpleName}"

    cpus = 1
    memory = 4.GB
    time = '1h'

    input:
    path(vcf)

    output:
    tuple path(vcf, includeInputs: true), path("${vcf}.tbi")

    script:
    """
    tabix ${vcf} 
    """

}

process MERGE {

    tag "${name}"

    cpus = 1
    memory = 4.GB
    time = '1h'

    input:
    path(vcf)
    path(tbi)
    val(name)

    output:
    path("${name}.vcf.bgz")

    script:
    """
    bcftools merge --output ${name}.vcf.bgz ${vcf}
    """
}

process MT {
    tag "${vcf.simpleName}"
    label 'hail'

    publishDir 'mt'

    cpus = 8
    memory = 16.GB
    time = '1h'

    input:
    path(vcf)

    output:
    path("${vcf.simpleName}.mt", type: 'dir')

    script:
    """
    #!python3

    import pandas as pd
    import hail as hl
    hl.init(local='local[${task.cpus}]')

    mt = hl.import_vcf("${vcf}",  reference_genome="GRCh38", array_elements_required=False)
    mt.write("${vcf.simpleName}.mt", overwrite=True)

    hl.stop()
    """
}