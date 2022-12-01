library(dplyr)
library(readr)

# this script removes the empty FASTA files from the blast results and parses the blast results to an updated config file
empty <- read_tsv(snakemake@input[['empty']], col_names = "query") %>%
  mutate(query = gsub("./", "", query),
         query = gsub(".fasta", "", query))
blast <- read_tsv(snakemake@input[['blast']],
                  col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
                  "qstart", "qend", "sstart", "send", "evalue", "bitscore"))

blast_empty_rm <- blast %>%
  filter(!sseqid %in% empty$query)

# late comers -- their seq download failed multiple times, then succeeded later.
# remove for now:
late_empty <- c("XP_014472329.1", 'KAG5370130.1', "XP_043497542.1", "XP_003041550.1",
                "OTA01579.1", "KAI5308474.1", "KAF6718986.1", "XP_027350160.1",
                "GBM10315.1", "XP_013320779.1", "KAI6854447.1", "AID23538.1")
blast_empty_rm <- blast_empty_rm %>%
  filter(!sseqid %in% late_empty)

# 47637 blast results after removing the empty FASTA sequences
length(unique(blast_empty_rm$sseqid))

# parse the results to a config file, which needs to be in the following format:
# query_protein:
#   - P63258_ACTG_BOVIN
#   - P51775_ACT_GIAIN
# features:
#   - atp_binding
#   - longitudinal_actin_contact
#   - lateral_actin_contact
features <- c("features:", "  - atp_binding", "  - lateral_actin_contact", "  - longitudinal_actin_contact")
config_string <- paste0("  - ", blast_empty_rm$sseqid)
config_string <- c("query_protein:", config_string, features)
write.table(config_string, file = "snakemake_config_blast.yml", 
            quote = F, row.names = F, col.names = F)
# I then back added the headers (query_protein, features) and the feature list.

# write the filtered blast results to a table
write_tsv(blast_empty_rm, snakemake@output[['out']])
