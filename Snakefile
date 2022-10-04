
configfile: "config.yml"

rule all:
    input:
        expand("{query_protein}.output.pdf", query_protein=config["query_protein"])

##############################################
##
##############################################

#####################################################
## Predicting Actin polymerization & ATPase activity
#####################################################

## This section uses the human beta-actin sequence as a reference
 
rule mafft_align_query_protein_to_human_beta_actin:
    input:
        human_beta_actin = 
        query_fasta = 
    output:
        mafft_map = 
        alignment=
    conda: "envs/mafft.yml"
    benchmark: "benchmarks/mafft_human_beta_actin_{query_protein}.txt"
    shell:'''
    mafft --auto --mapout --keeplength --add {input.query_fasta} {input.human_beta_actin} > {output.alignment}
    '''

rule calculate_shared_feature_residues:
    input: 
        feature_df = "inputs/protein_features/{features}.csv",
        mafft_map = ""
    output:
        df=,
        df_summary=
    conda: "envs/tidyverse.yml"
    benchmark: "benchmarks/calculate_shared_feature_residues_{query_protein}_{features}.txt"
    script: "snakemake/snakemake_calculate_shared_feature_residues.R"

##############################################
##
##############################################
