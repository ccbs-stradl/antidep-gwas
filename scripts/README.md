# Scripts for processing inputs or outputs of the nextflow scripts

- `make-bed_hg19.sh`: Creates separate bim, bed, fam files for each ancestry from the UK Biobank pgen files (hg19). **Check what nextflow pipeline output from this is used in...**

- `make-pgen_hg19.sh`: Use UKBiobank genetic data for creating LD references as these are in the hg19 format. Creates one pfile for all ancestries. With an extra column in the .psam file indicating the ancestry for that ID. Used in `fine_mapping.nf` and `genes.nf`.

- `fine_mapping_plots.R`: Formats and plots results from SuSiEx fine mapping. Used in `fine_mapping.nf`.

- `fixed.*`: Antidepressant exposure GWAS fixed-effects meta-analysis. Plots Manhattans and returns number of signifcant loci per genetic ancestry.

- `metasets.*`: Creates files with versions of each cohort used in each meta-analysis version eg. 2408 = meta-analysis run in August 2024, 2501 = meta-analysis run in Jan 2025.

- `multi.*`: Returns Manhattan plots and ancestry PCA plots from MR-MEGA. Also searches GWAS catalogue.

- `meta.pynb`: Meta-analysis using HAIL.
