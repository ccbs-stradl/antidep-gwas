# Antidepressant exposure GWAS: UK Biobank

Genome regression and association of antidepressant exposure phenotypes in [UK Biobank](https://www.ukbiobank.ac.uk).

## Setup

These scripts are primarily for documentation/reference purposes and file path setup is not fully described. Derivation of ATC codes requires [NHSBSA dm+d](https://isd.digital.nhs.uk/trud/user/guest/group/0/pack/6/subpack/24/releases) dictionary files downloaded from [NHS TRUD](https://isd.digital.nhs.uk/trud/). Phenotype derivation assumed local extracted copies of UKB data and is outdated.

## Scripts

- `antidepressants.R`: Create AD exposure phenotype from record-level health dataset stored in an [SQLite database](https://github.com/ccbs-stradl/ukb_healthoutcomes_db) (outdated)

## Workflows

- `ukb-top-regenie.nf`: Coordinate running of regenie between local cluster and RAP (outdated).