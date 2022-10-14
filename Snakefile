configfile: "snakemake_config.yml"

rule all:
    input:
        # PREDICTING FROM FEATURE RESIDUES
        expand("outputs/shared_feature_residues/1_shared_residue_information/{query_protein}-{features}.tsv", query_protein = config["query_protein"], features = config["features"]),
        expand("outputs/shared_feature_residues/2_shared_residue_summaries/{query_protein}-{features}.tsv", query_protein = config["query_protein"], features = config["features"]),
        # PREDICTING FROM HMM MODELS
        expand("outputs/hmm/hmmscan/{query_protein}-PF00022-hmmscan.out", query_protein = config["query_protein"])

#####################################################
## Predicting Actin polymerization & ATPase activity
#####################################################

# This section uses the human beta-actin sequence as a reference and hand-curated residue annotations for that protein.
# It determines the fraction of query residues for a given feature type that are shared with human beta-actin. 
 
rule mafft_align_query_protein_to_human_beta_actin:
    '''
    Aligns a query protein to the sequence for human beta actin. 
    The --auto flag allows the alignment algorith to be automatically selected; 
    in the context of two proteins, the most senstivie algorithm with usually be selected.
    Uses the --keeplength and --add parameters to enfore the reference coordinates to the alignment.
    This ensures that the curated residue annotations for the reference sequence can be transferred to the query protein.
    The --mapout flag outputs a textfile that maps the query protein into the reference sequence coordinates, 
    reporting the amino acid residue in the query, the query coordinate, and the reference coordinate.
    '''
    input:
        human_beta_actin = "inputs/P60709_ACTB_HUMAN.fasta",
        query_fasta = "query_proteins/{query_protein}.fasta"
    output:
        mafft_map = "query_proteins/{query_protein}.fasta.map",
        alignment= "outputs/shared_feature_residues/0_mafft/{query_protein}_vs_P60709_ACTB_HUMAN.fasta"
    conda: "envs/mafft.yml"
    benchmark: "benchmarks/mafft_human_beta_actin_{query_protein}.txt"
    shell:'''
    mafft --auto --mapout --keeplength --add {input.query_fasta} {input.human_beta_actin} > {output.alignment}
    '''

rule calculate_shared_feature_residues:
    '''
    Joins information in curated feature csv files to the mapout from mafft.
    Uses this information to calculcate the fraction of shared residues for each feature.
    '''
    input: 
        feature_csv = "inputs/protein_features/{features}.csv",
        mafft_map = "query_proteins/{query_protein}.fasta.map",
    output:
        tsv="outputs/shared_feature_residues/1_shared_residue_information/{query_protein}-{features}.tsv",
        tsv_summary="outputs/shared_feature_residues/2_shared_residue_summaries/{query_protein}-{features}.tsv"
    conda: "envs/tidyverse.yml"
    benchmark: "benchmarks/calculate_shared_feature_residues_{query_protein}_{features}.txt"
    script: "snakemake/snakemake_calculate_shared_feature_residues.R"

#####################################################
## Predict actin homology using hidden markov models
#####################################################

# This section uses a hidden markov model (HMM) built from PFAM actin alignments to score query protein queries for actin homology.
# The --cut_ga parameter passed to hmmscan filters homology predictions with low scores based on a profile-specific threshold.
# This translates the continuous measure of homology (E-value, score) into a binary classification: 
# only queries with strong homology to the hmm profile will pass the filtering threshold and be returned in the results.
# This section is built on PF00022, the PFAM profile for actin.
# For more information, see https://www.ebi.ac.uk/interpro/entry/pfam/PF00022/

rule download_pfam:
    output: "inputs/pfam/PF00022.hmm"
    benchmark: "benchmarks/hmm/PF00022-download-hmm.txt"
    shell:'''
    curl -JL --max-time 600 https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/PF00022?annotation=hmm | gunzip > {output}
    '''

rule hmmpress:
    input: "inputs/pfam/PF00022.hmm"
    output: "inputs/pfam/PF00022.hmm.h3i"
    conda: "envs/hmmer.yml"
    benchmark: "benchmarks/hmm/PF00022-hmmpress.txt"
    shell:'''
    hmmpress {input}
    '''

rule hmmscan:
    input:
        hmm = "inputs/pfam/PF00022.hmm",
        hmmpress = "inputs/pfam/PF00022.hmm.h3i",
        query_protein = "query_proteins/{query_protein}.fasta"
    output:
        out = "outputs/hmm/hmmscan/{query_protein}-PF00022-hmmscan.out",
        tbl = "outputs/hmm/hmmscan/{query_protein}-PF00022-hmmscan-tbl.out",
        dom = "outputs/hmm/hmmscan/{query_protein}-PF00022-hmmscan-dom.out"
    conda: "envs/hmmer.yml"
    benchmark: "benchmarks/hmm/hmmscan_{query_protein}.txt"
    shell:'''
    hmmscan --cut_ga -o {output.out} --tblout {output.tbl} --domtblout {output.dom} {input.hmm} {input.query_protein} 
    '''
