#!/bin/sh

sumstats=$1
output=$2

    #  1  CHROM
    #  2  POS
    #  3  ID
    #  4  A1
    #  5  AX
    #  6  CASE_ALLELE_CT
    #  7  CTRL_ALLELE_CT
    #  8  A1_CASE_FREQ
    #  9  A1_CTRL_FREQ
    # 10  MACH_R2
    # 11  FIRTH?
    # 12  TEST
    # 13  OBS_CT
    # 14  OR
    # 15  LOG(OR)_SE
    # 16  Z_STAT
    # 17  P
    # 18  N

gunzip -c $sumstats | awk 'OFS = "\t" {if(NR == 1) {print "chr", "pos", "ea", "oa", "beta", "se", "pval", "ncase", "ncontrol", "snp", "eaf", "imp_info", "eaf_case", "eaf_control"} else {print "chr"$1, $2, $4, $5, log($14), $15, $17, $6/2, $7/2, $3, $9, $10, $8, $9}}' > $output

# output columns
# chr
# pos
# ea
# oa
# beta
# se
# pval
# ncase
# ncontrol
# snp
# eaf
# imp_info
# eaf_case
# eaf_control