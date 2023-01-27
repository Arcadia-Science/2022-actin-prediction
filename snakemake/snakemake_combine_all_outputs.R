library(dplyr)
library(readr)

# read in average percent identity results
avg_pid <- read_tsv(snakemake@input[['avg_pid']]) %>%
  distinct() %>%
  # some tools inherit the seq name from the first line of the fasta,
  # and the accession needs to be edited to match the filename/genbank accession
  mutate(prot2 = gsub("sp\\|", "", prot2),
         prot2 = gsub("pdb\\|", "", prot2),
         prot2 = gsub("prf\\|\\|", "", prot2),
         prot2 = gsub("pir\\|\\|", "", prot2),
         prot2 = gsub("\\.1\\|.*",  "\\.1", prot2),
         prot2 = gsub("\\.2\\|.*",  "\\.2", prot2),
         prot2 = gsub("\\.3\\|.*",  "\\.3", prot2),
         prot2 = gsub("\\|", "_", prot2)) %>%
  filter(prot2 %in% blast_results$sseqid) %>%
  filter(!prot2 %in% c("NP_187818.1", "BAN39743.1", "TMW59056.1")) # filter three results where we detected inaccuracies
  
# read in the shared feature residue results (lateral contact, longitudinal contats, and ATP-binding residues)
features <- read_tsv(snakemake@input[['all_features']]) %>%
  distinct() %>%
  filter(protein %in% blast_results$sseqid)

# read in the hmm results
hmm <- read_tsv(snakemake@input[['all_hmm']]) %>%
  distinct() %>%
  mutate(query_name = gsub("sp\\|", "", query_name),
         query_name = gsub("pdb\\|", "", query_name),
         query_name = gsub("prf\\|\\|", "", query_name),
         query_name = gsub("pir\\|\\|", "", query_name),
         query_name = gsub("\\.1\\|.*",  "\\.1", query_name),
         query_name = gsub("\\.2\\|.*",  "\\.2", query_name),
         query_name = gsub("\\.3\\|.*",  "\\.3", query_name),
         query_name = gsub("\\|", "_", query_name))

# read in the foldseek results. The foldseek results also contain the uniprot accessions, protein names, and organism information
fsk <- read_tsv(snakemake@input[['all_fsk']]) %>%
  distinct() %>%
  mutate(evalue_transform = -1*log10(evalue))

# join everything together!
join1 <- full_join(features, avg_pid, by = c("protein" = "prot2"))
join2 <- full_join(join1, fsk, by = c("protein" = "genbank"))
all <- full_join(join2, hmm, by = c("protein" = "query_name")) %>%
  distinct()

write_tsv(all, snakemake@output[['all_outputs']])
