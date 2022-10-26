# Predicting whether a new protein fits into the conserved actin family using multiple computational criteria

## Purpose of this repository 

This repository uses a variety of computational approaches to predict whether a new protein is likely to be a functional actin sequence or not.

**Included approaches**

1. Fraction of specific functional residues shared between a query protein and an annotated actin protein.
2. Score query protein homology to PFAM actins with a hidden markov model.
3. Average percent identity between a query protein and a set of known actin proteins using pairwise protein alignments.

**Future approaches**

2. Structure comparisons (TBD)

## Getting started with this repository

This repository uses snakemake to run the pipeline and conda to manage software environments and installations.
You can find operating system-specific instructions for installing miniconda [here](https://docs.conda.io/en/latest/miniconda.html).
After installing conda and [mamba](https://mamba.readthedocs.io/en/latest/), run the following command to create the pipeline run environment.

```
mamba env create -n actin --file environment.yml
conda activate actin
```

To start the pipeline, run:
```
snakemake --use-conda -j 1
```

## Description of files and running your own data

The `Snakefile` orchestrates the pipeline that runs the computational approaches for actin prediction.
The pipeline is controlled by a config file, the path to which occurs on the first line of the `Snakefile`: `configfile = "snakemake_config.yml"`.
The config file provides two pieces of information: the functional predictions that should be executed (`atp_binding`, `longitudinal_actin_contact`, `lateral_actin_contact`) and the root name of the query protein sequences that the pipeline should be run on (`P63258_ACTG_BOVIN`).

To run your own protein of interest through the pipeline, copy the `snakemake_config.yml`, rename it, and replace the root name of the query protein sequence with the root name of your query protein sequence.
This should be everything that occurs before `*.fasta` in your protein sequence file name.
Put your protein sequence in the `query_proteins` directory.
Lastly, replace the path in the first line of the `Snakefile` with the path to your new config file and follow the "Getting started" instructions above.

Two config files are included in this repository: `snakemake_config.yml` and `snakemake_config_blast.yml`. 
`snakemake_config.yml` specifies two query protein sequences, and we have included these examples in this repository.
If you'd like to test the pipeline, you can run it with this config file by writing `configfile = "snakemake_config.yml"` on the first line of the `Snakefile`.
The `snakemake_config_blast.yml` specifies 50k additional protein query sequences. 
We obtained these sequences by running the `blast.snakefile` to download them.

## Developer instructions

### Snakemake

Snakemake manages rule-specific environments via the `conda` directive and using environment files in the `envs/` directory.
Snakemake itself is installed in the main conda environment `actin` specified in the `environment.yml` file.

Further details TBD.

### R code
To run unit tests on the R code, install the `testthat` package and libraries that are imported by the R functions.
```
mamba env create -n testthat --file testthat.yml
conda activate testthat
```

After installation, run the test script.
```
Rscript testthat.R
```
