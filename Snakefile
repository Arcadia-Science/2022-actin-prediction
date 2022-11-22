configfile: "snakemake_config_blast.yml"

rule all:
    input:
        # ESTIMATING PAIRWISE IDENTITY
        expand("outputs/mean_pid/{query_protein}_pid.tsv", query_protein = config["query_protein"]),
        # PREDICTING FROM FEATURE RESIDUES
        expand("outputs/shared_feature_residues/1_shared_residue_information/{query_protein}-{features}.tsv", query_protein = config["query_protein"], features = config["features"]),
        expand("outputs/shared_feature_residues/2_shared_residue_summaries/{query_protein}-{features}.tsv", query_protein = config["query_protein"], features = config["features"]),
        # PREDICTING FROM HMM MODELS
        expand("outputs/hmm/hmmscan/{query_protein}-PF00022-hmmscan.out", query_protein = config["query_protein"])

#####################################################
## Estimating average pairwise identity between a 
## query protein and known actin protein sequences
#####################################################

# This section calculates pairwise identity between a query protein and a curated set of actin proteins that are confidentally actin protein sequences.

rule create_fasta_file_with_confident_actins_and_query_proteins:
    input:
        confident_fasta="inputs/confident_actins.fasta",
        query_fasta = "query_proteins/{query_protein}.fasta"
    output: "outputs/mean_pid/{query_protein}_combined.fasta"
    benchmark: "benchmarks/combine_with_confident_{query_protein}.txt"
    shell:'''
    cat {input.confident_fasta} {input.query_fasta} > {output}
    ''' 

rule mafft_multiple_sequence_align_query_protein_and_confident_actins:
    '''
    This rule generates a multiple sequence alignment between confident actins and a query protein sequence.
    It uses the program mafft-linsi, an alias for an accurate option (L-INS-i) for an alignment of up to ∼200 sequences × ∼2,000 sites.
    '''
    input: "outputs/mean_pid/{query_protein}_combined.fasta"
    output: "outputs/mean_pid/{query_protein}_msa.fasta"
    conda: "envs/mafft.yml"
    benchmark: "benchmarks/msa_with_confident_actins_{query_protein}.txt"
    shell:'''
    mafft-linsi {input} > {output}
    '''

rule calculate_pairwise_identity:
    input: 
        msa = "outputs/mean_pid/{query_protein}_msa.fasta"
    output: 
        tsv = "outputs/mean_pid/{query_protein}_pid.tsv",
        pdf = "outputs/mean_pid/{query_protein}_pid.pdf",
    conda: "envs/tidybio3d.yml"
    benchmark: "benchmarks/calculate_pid_{query_protein}.txt"
    script: "snakemake/snakemake_calculate_pairwise_identity.R"

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

#########################################################
## Compare structure of query protein against known actin
#########################################################

checkpoint chunk_genbank_accessions_for_uniprot_id_conversion:
    '''
    checkpoints instruct snakemake to re-evaluate the DAG mid-workflow.
    This allows us to create a new wildcard from information generated by the workflow itself. 
    In this case, we'll dynamically determine how many "chunks" we need from our input accessions to limit each uniprot API query to 15k or fewer protein accessions.
    '''
    input: fastas=expand("query_proteins/{query_protein}.fasta", query_protein = config["query_protein"])
    output: outdir=directory("outputs/foldseek/chunked_accessions/")
    conda: "envs/tidyverse.yml"
    benchmark: "benchmarks/foldseek/chunk_query_proteins.txt"
    script: "snakemake/snakemake_chunk_genbank_accessions_for_uniprot_id_conversion.R" 

rule make_query_form: 
    """
    This rule formats the genbank accessions into a uniprot API query
    """
    input: "outputs/foldseek/chunked_accessions/genbank_accessions_chunk{grp}.txt"
    output:
        query_form="outputs/foldseek/uniprot_accessions/query_form{grp}.txt",
    benchmark: "benchmarks/foldseek/convert_genbank_to_uniprot{grp}.txt"
    shell:"""
    # store the genbank ids as a variable:
    ids=$(cat {input} | tr "\n" ",") # replace newlines with a comma

    # write the curl forms to a file - this is needed otherwise specifying the IDs as a variable will fail:
    echo "from=EMBL-GenBank-DDBJ_CDS&to=UniProtKB&ids=$ids" > {output.query_form}
    """

rule convert_genbank_protein_accession_to_uniprot_accessions:
    """
    This rule converts a genbank accessions to uniprot using the api query created in make_query_form.
    It submits the query to the api, checks the status of the query, and then downloads the results when the query is done running
    This will be lossy -- not every genbank accession will have a uniprot accession.
    This could be broken up into 2 rules, but I think I like it better as one?
    The first would make the query and the second would create the job id, check the status, and download results
    I like it as one so that the query submission is coupled with the download. 
    """
    input: query_form="outputs/foldseek/uniprot_accessions/query_form{grp}.txt",
    output:
        query="outputs/foldseek/uniprot_accessions/query{grp}.txt",
        query_status="outputs/foldseek/uniprot_accessions/query_status{grp}.txt",
        results="outputs/foldseek/uniprot_accessions/results{grp}.txt"
    conda: "envs/curl.yml"
    benchmark: "benchmarks/foldseek/convert_genbank_to_uniprot{grp}.txt"
    shell:"""
    # submit the curl "form" to uniprot to start the query
    curl -d @{input.query_form} https://rest.uniprot.org/idmapping/run > {output.query}

    # extract the job ID from this query
    jobID=$(cat {output.query} | sed "s/.*://g" | sed "s/}}//g" | sed 's/"//g')

    waiting=1
    finished='{{"jobStatus":"FINISHED"}}'
    while [ $waiting == 1 ]
        do
        # wait a few seconds
        sleep 10

        # download the query status
        curl -i "https://rest.uniprot.org/idmapping/status/$jobID" | tail -n1 > {output.query_status}
        qstatus="$(cat {output.query_status})"
        echo $qstatus

        # check whether the uniprot query is still running. If so keep going, otherwise download the results.
        if [[ "$qstatus" != "$finished" ]]
        then
            waiting=1 # Still waiting
        else
            waiting=0
            echo "downloading results"
            #curl -s "https://rest.uniprot.org/idmapping/uniprotkb/results/stream/$jobID?format=tsv" > {output.results} # download the mapped ids.
            curl -JLo {output.results} https://rest.uniprot.org/idmapping/uniprotkb/results/stream/$jobID?fields=accession%2Cid%2Creviewed%2Cprotein_name%2Cgene_names%2Corganism_name%2Clength&format=tsv
        fi
    done
    """

def chunk_genbank_accessions_for_uniprot_id_conversion(wildcards):
    # expand checkpoint to get grp values, and place them in the final file name that uses that wildcard
    # checkpoint_output encodes the output dir from the checkpoint rule. 
    checkpoint_output = checkpoints.chunk_genbank_accessions_for_uniprot_id_conversion.get(**wildcards).output[0]    
    file_names = expand("outputs/foldseek/uniprot_accessions/results{grp}.txt",
                        grp = glob_wildcards(os.path.join(checkpoint_output, "genbank_accessions_chunk{grp}.txt")).grp)
    return file_names

rule combine_uniprot_id_conversions:
    input: results=chunk_genbank_accessions_for_uniprot_id_conversion,
    output: tsv = "outputs/foldseek/uniprot_accessions/results.tsv"
    conda: "envs/tidyverse.yml"
    benchmark: "benchmarks/foldseek/convert_genbank_to_uniprot_combine_results.txt"
    script: "snakemake/snakemake_combine_uniprot_id_conversions.R"

checkpoint create_dummy_files_for_uniprot_accession_wildcard:
    input: tsv = "outputs/foldseek/uniprot_accessions/results.tsv"
    output: outdir = directory("outputs/foldseek/uniprot_accessions_wc/")
    script: "snakemake/snakemake_create_dummy_files_for_uniprot_accession_wildcard.R"

rule download_alphafold_pdb_files_for_uniprot_accessions:
    input: "outputs/foldseek/uniprot_accessions_wc/{uniprot_acc}.txt"
    output: "outputs/foldseek/uniprot_alphafold_pdb_structures/AF-{uniprot_acc}-F1-model_V4.pdb"
    conda: "envs/curl.yml"
    benchmark: "benchmarks/foldseek/download_alphafold_pdb_files/{uniprot_acc}.txt"
    shell:'''
    curl -JLo {output} https://alphafold.ebi.ac.uk/files/AF-{wildcards.uniprot_acc}-F1-model_v4.pdb
    '''

def create_dummy_files_for_uniprot_accession_wildcard(wildcards):
    # expand checkpoint to get uniprot acc values, and place them in the final file name that uses that wildcard
    # checkpoint_output encodes the output dir from the checkpoint rule. 
    checkpoint_output = checkpoints.create_dummy_files_for_uniprot_accession_wildcard.get(**wildcards).output[0]    
    file_names = expand("outputs/foldseek/foldseek/{uniprot_acc}_vs_1j6z.tsv",
                        uniprot_acc = glob_wildcards(os.path.join(checkpoint_output, "{uniprot_acc}.txt")).uniprot_acc)
    return file_names

rule run_foldseek:
    """
    The reference sequence is rabbit actin and was downloaded from this web page: https://www.rcsb.org/structure/1J6Z
    """
    input:
        reference = "inputs/pdb/1j6z.pdb",
        query = "outputs/foldseek/uniprot_alphafold_pdb_structures/AF-{uniprot_acc}-F1-model_V4.pdb"
    output: "outputs/foldseek/foldseek/{uniprot_acc}_vs_1j6z.tsv"
    params: refdir = "inputs/pdb/"
    conda: "envs/foldseek.yml"
    shell:'''
    foldseek easy-search {input.query} {params.refdir} {output} tmp_foldseek_folder 
    '''

rule tmp:
    input: create_dummy_files_for_uniprot_accession_wildcard
