//nextflow.enable.dsl=2

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
	
	// Meta-analysis
	MA_CH = MA(MT_CH)
	
	ANALYSES_CH = MA_CH
		.map { it -> it[1] }
		.splitCsv(header: true)
		.map { it -> it.analysis }
	
	// export
	TSV_CH = EXPORT(MA_CH, ANALYSES_CH)
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

process MA {
	tag "${mt.simpleName}"
	label 'hail'
	
	cpus = 8
	memory = 16.GB
	time = '1h'
	
	input:
	path(mt)
	
	output:
	tuple path("${mt.simpleName}-ma.mt", type: 'dir'), path("analyses.csv")
	
	script:
	"""
	#!python3
	
	import pandas as pd
	import hail as hl
	hl.init(local='local[${task.cpus}]')
	
	# read sumstats
	gw = hl.read_matrix_table("${mt}")
	
	# annotate datasets with analysis grouping (phenotype + cluster)
	samples = gw.s.collect()
	meta_df = pd.DataFrame(
		{
		"dataset": samples,
		"analysis": ["-".join(s.split('-')[1:3]) for s in samples]
		}
	)
	meta_t = hl.Table.from_pandas(meta_df, key="dataset")
	
	gw = gw.annotate_cols(analysis=meta_t[gw.s].analysis)
	
	# meta-analysis weights
	
	# inverse variance weight for each estimate
	# sample size weight for other data
	mw = gw.annotate_entries(
		ivw=1 / gw.SE ** 2,
		ncw=1 / gw.NC
	)
	
	# effect weighted by inverse variance
	mw = mw.annotate_entries(
		ES_w=mw.ES * mw.ivw,
		AF_w=mw.AF * mw.ncw,
		SI_w=mw.SI * mw.ncw
	)
	
	# meta-analysis aggregation
	ma = mw.group_cols_by(mw.analysis).aggregate(
		ES=(hl.agg.array_sum(mw.ES_w) / hl.agg.array_sum(mw.ivw)).first(),
		SE=hl.map(lambda ivw_sum: hl.sqrt(1 / ivw_sum), hl.agg.array_sum(mw.ivw)).first(),
		NC=hl.agg.array_sum(mw.NC).first(),
		SS=hl.agg.array_sum(mw.SS).first(),
		AF=(hl.agg.array_sum(mw.AF_w) / hl.agg.array_sum(mw.ncw)).first(),
		SI=(hl.agg.array_sum(mw.SI_w) / hl.agg.array_sum(mw.ncw)).first()
	)
	
	ma = ma.annotate_entries(
		PP=hl.pchisqtail(ma.ES ** 2 / ma.SE ** 2, df=1)
	)
	
	ma.write("${mt.simpleName}-ma.mt")
	
	# analyses list
	analyses = pd.DataFrame({"analysis": ma.analysis.collect()})
	analyses.to_csv("analyses.csv", index=False)
	
	hl.stop()
	"""
}

process EXPORT {
	tag "${analysis}"
	label 'hail'
	
	cpus = 8
	memory = 16.GB
	time = '1h'
	
	input:
	tuple path(ma), path(csv)
	each analysis
	
	output:
	path("hail-${analysis}.tsv.bgz")
	
	script:
	"""
	#!python3
	
	import pandas as pd
	import hail as hl
	hl.init(local='local[${task.cpus}]')
	
	# read sumstats
	ma = hl.read_matrix_table("${ma}")
	
	# write out flat table for one analysis
	ma.filter_cols(ma.analysis == "${analysis}").make_table().flatten().export("hail-${analysis}.tsv.bgz")
	
	"""
}