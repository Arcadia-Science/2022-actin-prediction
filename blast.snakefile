rule all:
    input: "outputs/blast/blastp_results.out"

rule blast:
    input: "inputs/P60709_ACTB_HUMAN.fasta"
    output: "outputs/blast/blastp_results.out"
    conda: "envs/blast.yml"
    shell:'''
    blastp -db nr -query {input} -out {output} -remote -max_target_seqs 50000 -outfmt 6
    '''
