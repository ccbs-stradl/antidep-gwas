# Scripts for processing inputs or outputs of the nextflow scripts

- `make-bed_hg19.sh`: Creates separate bim, bed, fam files for each ancestry from the UK Biobank pgen files (hg19). **Check what nextflow pipeline output from this is used in...**

- `make-pgen_hg19.sh`: Use UKBiobank genetic data for creating LD references as these are in the hg19 format. Creates one pfile for all ancestries. With an extra column in the .psam file indicating the ancestry for that ID. Used in `fine_mapping.nf` and `genes.nf`.

- `fine_mapping_plots.R`: Formats and plots results from SuSiEx fine mapping. Used in `fine_mapping.nf`.