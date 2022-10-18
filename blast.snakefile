import pandas as pd
import os

rule all:
    input: "qpa_done.txt"

rule blast:
    input: "inputs/P60709_ACTB_HUMAN.fasta"
    output: "outputs/blast/blastp_results.out"
    conda: "envs/blast.yml"
    shell:'''
    blastp -db nr -query {input} -out {output} -remote -max_target_seqs 50000 -outfmt 6
    '''

checkpoint download_query_protein_fastas:
    input: "outputs/blast/blastp_results.out"
    output: directory("outputs/blast/query_protein_fastas/")
    run:
        # check if output director exists, and if not, create it
        outdir = output[0]
        check_outdir_exists = os.path.isdir(outdir)
        if not check_outdir_exists:
            os.makedirs(outdir)

        # read in blast results, iterate through the matched protein accessions, and download as fastas
        blastp_results = pd.read_csv(str(input[0]), header = None, sep = "\t")
        query_protein_accessions = blastp_results[1].unique().tolist()
        for query_protein_accession in query_protein_accessions:
            output_fasta_file = query_protein_accession + ".fasta"
            output_fasta_path = os.path.join("outputs", "blast", "query_protein_fastas", output_fasta_file) 
            shell('esearch -db protein -query "{query_protein_accession}" | efetch -format fasta > {output_fasta_path}')

def checkpoint_download_query_protein_fastas(wildcards):
    # Expand checkpoint to get fasta names, which will be used as query proteins for the main Snakefile
    # checkpoint_output encodes the output dir from the checkpoint rule.
    checkpoint_output = checkpoints.download_query_protein_fastas.get(**wildcards).output[0]    
    file_names = expand("outputs/blast/query_protein_fastas/{query_protein_accession}.fasta",
                        query_protein_accession = glob_wildcards(os.path.join(checkpoint_output, "{query_protein_accession}.fasta")).query_protein_accession)
    return file_names


rule collapse_tmp:
    input: checkpoint_download_query_protein_fastas
    output: touch("qpa_done.txt")
