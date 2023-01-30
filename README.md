# Predicting whether a new protein fits into the conserved actin family using multiple computational criteria

## Purpose of this repository 

This repository uses a variety of computational approaches to predict whether a new protein is likely to be a functional actin sequence or not.

**Included approaches**

1. Fraction of specific functional residues shared between a query protein and an annotated actin protein.
2. Score query protein homology to PFAM actins with a hidden markov model.
3. Average percent identity between a query protein and a set of known actin proteins using pairwise protein alignments.
4. Structure comparisons between a query protein (alphafold-predicted structures) and rabbit beta actin using foldseek.

## How to run your own FASTA file through the pipeline using Binder

TODO: update binder branch to `binder`
[![Binder](https://aws-uswest2-binder.pangeo.io/badge_logo.svg)](https://aws-uswest2-binder.pangeo.io/v2/gh/Arcadia-Science/2022-actin-prediction/ter/binderize)

We built a [Binder](https://mybinder.org/) to enable others to run the actin pipeline on there protein FASTA sequences.
If you're interested in running the pipeline on your sequence(s) of interest, make sure they meet the following guidelines:

1. The protein sequence needs to be in amino acid FASTA format.
2. The sequence should represent a suspected actin, actin-like, actin-related, or other sequence that is similar to actin.
3. The sequence must be available on GenBank. This is because the pipeline uses pre-calculated alphafold PDB structures based on the UniProt database. The pipeline provides a conversion between the GenBank identifier and the UniProt identifier, but the input file name must be the GenBank identifer (see below).

You need to do three things start the pipeline on your own FASTA sequence.
We provide an overview of these steps first, and then provide more details below.

1. Launch Binder ([![Binder](https://aws-uswest2-binder.pangeo.io/badge_logo.svg)](https://aws-uswest2-binder.pangeo.io/v2/gh/Arcadia-Science/2022-actin-prediction/ter/binderize))
2. Upload your FASTA sequence to the Binder and edit the config file to tell the pipeline the name of your FASTA file.
3. Start the snakemake pipeline.

### 1. Launch Binder

The first step is to launch the Binder instance.
Binder is a service that turns a Git repo into a collection of interactive notebooks that can be executed on a temporary cloud machine.
Binder is currently free, so there is no cost to use to try it out!
We built the Binder for this repository on [Pangeo Binder](https://pangeo-binder.readthedocs.io/en/prod/).
This BinderHub provides slightly more powerful compute, but [to combat users who were using the service illicitly](https://github.com/pangeo-data/pangeo-binder/issues/195), you have to login via GitHub and can only have one instance running at a time. 
GitHub accounts are also free, so you can create one if you don't already have one.

<details>
  <summary>More information on binder and what happens when you click the launch binder button.</summary>
Binder is a service that turns a Git repo into a collection of interactive notebooks. 
When a repository is configured to run as a binder, passing the GitHub repository URL to binder starts the binder-building process.
Binder first builds a docker image that contains all of the software installations specified by a special set of files in the GitHub repository.
A docker image is a set of instructions that are used to create a docker container.
A docker container is a runnable instance of a docker image -- it's an encapsulated computing environment that can be used to reproducibly install sets of software on diverse computers.
Armed with the docker container, binder launches an "instance" in the cloud (either on Google Cloud or AWS typically) on which it runs the docker container.
Binder does some additional work in the background -- if no software configuration files are provided in the GitHub repo, or if those contain a minimal set of software, binder will by default include JupyterHub in the docker.
When the cloud instance is launched, this is the screen you interact with.
You interact with the cloud instance in your browser.
Binders are ephemeral instances -- after a period of inactivity, the instance is automatically shut down, and any work you have done will be lost.
You're able to download files from your work before the instance is shut down if you do want to save anything.
</details>
 
### Upload your FASTA sequence to the Binder and edit the config file to tell the pipeline the name of your FASTA file

Once you have the Binder running, you need to upload your FASTA file to the Binder.
The file needs to be uploaded into the `query_proteins` folder.
Double click on this folder to enter into it.

![](https://i.imgur.com/0vnfjZi.png)

Then, use the upload arrow to get a dialogue prompt that points to files on your local computer to upload the FASTA sequence.



* if there is no foldseek output, those columns won't be in the output CSV file

Two config files are included in this repository: `snakemake_config.yml` and `snakemake_config_blast.yml`. 
`snakemake_config.yml` specifies two query protein sequences, and we have included these examples in this repository.
If you'd like to test the pipeline, you can run it with this config file by writing `configfile = "snakemake_config.yml"` on the first line of the `Snakefile`.

To run your own protein of interest through the pipeline, copy the `snakemake_config.yml`, rename it, and replace the root name of the query protein sequence with the root name of your query protein sequence.
This should be everything that occurs before `*.fasta` in your protein sequence file name.
Put your protein sequence in the `query_proteins` directory.
Lastly, replace the path in the first line of the `Snakefile` with the path to your new config file and follow the "Getting started" instructions above.

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

## Description of files & how the pipeline works

The `Snakefile` orchestrates the pipeline that runs the computational approaches for actin prediction.
The pipeline is controlled by a config file, the path to which occurs on the first line of the `Snakefile`: `configfile = "snakemake_config.yml"`.
The config file provides two pieces of information: the functional predictions that should be executed (`atp_binding`, `longitudinal_actin_contact`, `lateral_actin_contact`) and the root name of the query protein sequences that the pipeline should be run on (`P63258_ACTG_BOVIN`).

In the `main` branch of this repository, the `Snakefile` is run with the config file `snakemake_config_blast.yml`.
This config file specifies ~50k protein query sequences that we ran through the pipeline. 
We obtained these sequences by running the [`blast.snakefile`](./blast.snakefile) to download them.
Note that due to an error in the NCBI entrez API (`Invalid uid * at position= 0`), only 48,406 of the 50k sequences were downloadable, so only these were included in our analysis. 
The Snakefile in the `main` branch is optimized for the `snakemake_config_blast.yml` config file, including working with the file output by the BLAST run and producing plots that are appropriate for thousands of results.
A generalized `Snakefile` that will work with new actin FASTA files and doesn't rely on BLAST results is available in the [`binder`](https://github.com/Arcadia-Science/2022-actin-prediction/tree/binder) branch.
See the [above section](#how-to-run-your-own-fasta-file-through-the-pipeline-using-binder) for instructions on how to run a new actin FASTA sequence using that Snakefile.

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
