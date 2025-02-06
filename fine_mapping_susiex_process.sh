CLUMP_RANGES_FILE=${finemapRegions}

tail -n +2 "${CLUMP_RANGES_FILE}" | while IFS=$'\t' read -r line; do
  # Get CHR, BP_START, BP_END by using awk to match column names indices
  # Really important we put these cols in the correct order in the R scripts making the "loci" object
  CHR=$(echo "$line" | awk '{print $1+0}')
  BP_START=$(echo "$line" | awk '{print $2+0}')
  BP_END=$(echo "$line" | awk '{print $3+0}')

  echo "Processing CHR: $CHR, BP_START: $BP_START, BP_END: $BP_END"

  SuSiEx \ 
    --sst_file=${maPaths} \
    --n_gwas=${neff} \
    --ref_file=${bfile} \
    --ld_file=$(echo "${ancestries}" | sed 's/\\([^,]*\\)/\\1.${chr}/g') \
    --out_dir=. \
    --out_name=SuSiEx.${ancestries}.output.cs95_${CHR}:${BP_START}:${BP_END} \
    --level=0.95 \
    --pval_thresh=1e-5 \
    --chr=$CHR \
    --bp=$BP_START,$BP_END \
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
    --keep-ambig=True |& tee SuSiEx.${ancestries}.output.cs95_${CHR}:${BP_START}:${BP_END}.log

  break

done 
