# GWAS Meta-analysis of Antidepressant Prescribing

This repository contains scripts, documentation and configuration files for conducting a genome-wide association study (GWAS) meta-analysis of antidepressant exposure.

## Repository Structure

- `docs/`: Contains documentation for running the pipelines
    - `workflow_documentation.md` step-by-step instructions for setting up the environment, downloading necessary files, and running the various Nextflow pipelines to reproduce the GWAS meta-analysis of antidepressant exposure.
    - `workflow_references.md` instructions to download required reference files for the workflows.
    - `workflow_requirements.md` indsutrctions for installing any requirements.
- `workflows/` : Contains nextflow scripts for running the pipelines.
- `scripts/`: Contains R markdown scripts (and their outputs), jupyter notebooks, bash shell scripts for data pre-processeing and formatting/plotting outputs from the nextflow scripts. See the README.md in this directory for more info on what each script does.
- `results/`: Contains outputs from the nextflow pipeline scripts. Sub-folders are made inside this directory for each nextflow pipeline. Most of these subfolders are not git commited, instead symlinks should be made to point to their contents (see `docs/workflow_documentation.md`)
- `config/`: Contains nextflow configuration files
- `env/`: Contains yaml files for setting up virtual conda environments.
- `manuscript/` : Contains scripts and outputs used directly for the manuscript, eg. creating tables and pasting sentences with sample size numbers that can be copied into the main text of the manuscript.
- `reference/` : Contains reference files needed in nextflow scripts, see `docs/workflow_references.md`


## Contact

For questions or collaboration inquiries, please contact Mark Adams (mark.adams@ed.ac.uk) and Amelia Edmondson-Stait (amelia.edmondson-stait@ed.ac.uk).



