# Predicting whether a new protein fits into the conserved actin family using multiple computational criteria

## Purpose of this repository 

This repository uses a variety of computational approaches to predict whether a new protein is likely to be a functional actin sequence or not.

**Included approaches**

1. Fraction of specific functional residues shared between a query protein and an annotated actin protein.

**Future approaches**

1. Average percent identity between a query protein and a set of known actin proteins using pairwise protein alignments.
2. Score with a hidden markov model using PFAM actin hmm. 
3. Structure comparisons (TBD)

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
