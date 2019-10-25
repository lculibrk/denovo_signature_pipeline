## Mutation Signature Pipeline

### Getting started
Clone

    git clone https://github.com/lculibrk/denovo_signature_pipeline/

Build a conda environment from the conda spec file and activate it:

    conda create -p ./signature_env --file conda_spec.txt
    conda activate signature_env/

Place data in .maf format under data/{project}/genome/cohorts/{cohort}/snv/{sampleid}.maf

Edit the config.yaml file with your specific parameters

Run Snakemake:

    snakemake
    
Alternatively use the launch.sh as a template (or as-is on the GSC servers) to run the pipeline on a cluster:

    nohup bash launch.sh 1000 Snakefile &
    
