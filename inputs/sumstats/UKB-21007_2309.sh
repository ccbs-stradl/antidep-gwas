#!/bin/sh

sumstats=$1
output=$2

# sumstats columns
#  1  CHROM
#  2  GENPOS
#  3  ID
#  4  ALLELE0
#  5  ALLELE1
#  6  A1FREQ
#  7  A1FREQ_CASES
#  8  A1FREQ_CONTROLS
#  9  INFO
# 10  N
# 11  N_CASES
# 12  N_CONTROLS
# 13  TEST
# 14  BETA
# 15  SE
# 16  CHISQ
# 17  LOG10P
# 18  EXTRA

gunzip -c $sumstats | awk 'OFS = "\t" {if(NR == 1) {print "chr", "pos", "ea", "oa", "beta", "se", "pval", "ncase", "ncontrol", "snp", "eaf", "imp_info", "eaf_case", "eaf_control", "neff"} else {print $1, $2, $5, $4, $14, $15, 10^(-$17), $11, $12, $3, $8, $9, $7, $8, 4 / (1/$11 + 1/$12)}}' > $output

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
# neff