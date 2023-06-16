# Predicting whether a new protein fits into the conserved actin family using multiple computational criteria

## Purpose of this repository 

This repository uses a variety of computational approaches to predict whether a new protein is likely to be a functional actin sequence or not.

**Included approaches**

1. Fraction of specific functional residues shared between a query protein and an annotated actin protein.
2. Score query protein homology to PFAM actins with a hidden markov model.
3. Average percent identity between a query protein and a set of known actin proteins using pairwise protein alignments.
4. Structure comparisons between a query protein (alphafold-predicted structures) and rabbit beta actin using foldseek.

## How to run your own FASTA file through the pipeline using Binder

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/Arcadia-Science/2022-actin-prediction/binder)

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
 
We built a [Binder](https://mybinder.org/) to enable others to run the actin pipeline on there protein FASTA sequences.
If you're interested in running the pipeline on your sequence(s) of interest, make sure they meet the guidelines below.
An overview is provided here with more detailed instructions below.

1. The protein sequence needs to be in amino acid FASTA format.
2. The sequence should represent a suspected actin, actin-like, actin-related, or other sequence that is similar to actin.
3. The sequence must be available on GenBank. This is because the pipeline uses pre-calculated alphafold PDB structures based on the UniProt database. The pipeline provides a conversion between the GenBank identifier and the UniProt identifier, but the input file name must be the GenBank identifer (see below).

You need to do three things to start the pipeline on your own FASTA sequence.
We provide an overview of these steps first, and then provide more details below.

1. Launch Binder [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/Arcadia-Science/2022-actin-prediction/binder)
2. Upload your FASTA sequence to the Binder and edit the config file to tell the pipeline the name of your FASTA file.
3. Start the snakemake pipeline.

### 1. Launch Binder

The first step is to launch the Binder instance.
Binder is a service that turns a Git repo into a collection of interactive notebooks that can be executed on a temporary cloud machine.
Binder is currently free, so there is no cost to using it!

### 2. Upload your FASTA sequence to the Binder and edit the config file to tell the pipeline the name of your FASTA file

Once you have the Binder running, you need to upload your FASTA file to the Binder.
The file needs to be uploaded into the `query_proteins` folder.
Double click on this folder to enter into it.

![](https://i.imgur.com/0vnfjZi.png)

Then, use the upload arrow to get a dialogue prompt that points to files on your local computer to upload the FASTA sequence.
Upload your file(s).

![](https://i.imgur.com/wUzgLFv.png)

Next, edit the config file so that it shows the GenBank accession number of your protein of interest.
Start by double clicking on the config file so that it opens in a text editor.

![](https://i.imgur.com/bzIHtI7.png)

Then, edit the file so that your accession number is listed under `query_protein:`.
The accession number needs to be the same as the file name prefix that you uploaded (see the test data in `query_proteins` and annotated in the config file `snakemake_config.yml` for an example).
Save the file and exit from the text editor.

### 3. Start the snakemake pipeline

The last step is to start the Snakemake pipeline.
From the Launcher tab, use the `Terminal` icon to open a new terminal.

![](https://i.imgur.com/SozBrHJ.png)

Then, enter the command `snakemake -j 2` and allow the pipeline to run.
If you uploaded a single protein, it should take about 1-2 minutes.
As the pipeline starts, you should see something like:

![](https://i.imgur.com/4R7x77J.png)

### Explanation of outputs

When the pipeline is finished, you'll have a new file named `outputs/all_outputs_summarized.tsv`. 
If you double click on it, it will look something like this:

![](https://i.imgur.com/Regpc8b.png)

Below we provide a description of the columns:

* `protein`: GenBank accession of protein.
* **Functional annotation**: functional annotation based on residue identity at residues involved in the polymerization and ATPase activity of actin. `lon` stands for contacts that are important for longitudinal polymerization, `lat` stands for contacts that are important for latitudinal polymerization, and `atp` stands for residues that are important for ATPase activity.
  * `*_feature_count`: the number of residues that were important for a specific function.
  * `*_num_matching`: the number of residues that matched the reference.
  * `*_fraction_matching`: the fraction of all possible matches within the category that matched to the reference.
  * `w_avg_contact`: weighted average for all categories `lon`, `lat`, and `atp`.
* **Alignment**: alignment annotation based on query alignment to a multiple sequence alignment of reference actin sequences.
  * `avg_fid`: average fraction identity.
  * `avg_pid`: average percent identity.
* **Structure**: structural annotation performed by comparing the alphafold-predicted query structure against the PDB structure of monomeric rabbit actin. If there is no Foldseek output, those columns won't be in the output CSV file. See the [Foldseek preprint](https://www.biorxiv.org/content/10.1101/2022.02.07.479398v1) for more information on interpretting these output columns.
  * `query`: UniProt ID
  * `target`: reference PDB structure used for structural comparison (rabbit actin).
  * `pident`: percent sequence identity between the query and target sequences.
  * `alnlen`: alignment length
  * `mismatch`: number of mismatches in the alignment.
  * `gapopen`: gaps in alignment.
  * `qstart`: start residue for query in alignment
  * `qend`: end residue for query in alignment
  * `tstart`: start residue for target in alignment
  * `tend`: end residue for target in alignment
  * `evalue`: main output of Foldseek that reflects the structural similarity of the query and target proteins.
  * `bits`: information measurement
  * `uniprot_name`: name of protein in UniProt
  * `Reviewed`: reviewed status of protein
  * `Protein names`: UniProt protein name
  * `Gene Names`: UniProt gene name
  * `Organism`: UniProt organism
  * `Length`: UniProt length of protein
  * `evalue_transform`: `-1*log10(evalue)`
* **Hidden markov model**: sequence similarity annotation performed using a hidden markov model of PFAM actin. See the [HMMER user guide](http://eddylab.org/software/hmmer3/3.1b2/Userguide.pdf) for more information about the output columns.
  * `domain_name`: PFAM name
  * `domain_accession`: PFAM accession
  * `query_accession`: query accession in PFAM (`NA`)
  * `sequence_evalue`: evalue of sequence compared to the PFAM model; statistical significance of the query to the PFAM model.
  * `sequence_score`: log-odds score of sequence compared to PFAM model. Higher scores are more similar to the HMM model.
  * `sequence_bias`: correction term for biased sequence composition.
  * `best_domain_evalue`: same as above, but for best-scoring domain in the sequence (some accessions have multiple domains).
  * `best_domain_score`: same as above
  * `best_domain_bis`: same as above
  * `domain_number_exp`: the posterior probability of alignment starts and ends (profile B and E state alignments) with respect to target sequence position. The sum of the posterior probabilities of alignment starts (B states) over the entire target sequence is the expected number of domains in the sequence.
  * `domain_number_reg`: number of discrete regions identified by this posterior decoding step. 
  * `domain_number_clu`: the number of regions that had to be subjected to stochastic traceback clustering.
  * `domain_number_ov`: for envelopes that were defined by stochastic traceback clustering, how many of them overlap other envelopes.
  * `domain_number_env`: total number of domain envelopes identified.
  * `domain_number_dom`: number of domains defined.
  * `domain_number_rep`: number of domains satisfying reporting thresholds.
  * `domain_number_inc`:  number of domains satisfying inclusion thresholds.
  * `description`: description of the target, as free text.

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
