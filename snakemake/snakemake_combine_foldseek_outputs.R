library(readr)
library(dplyr)
library(purrr)

fsk_files <- unlist(snakemake@input[['fsk']])

fsk <- fsk_files %>%
  purrr::map_dfr(readr::read_tsv, show_col_types = F, 
                 col_names = c("query", "target", "pident", "alnlen", "mismatch", 
                               "gapopen", "qstart", "qend", "tstart", "tend", "evalue", 
                               "bits"),
                 col_types = "ccddddddddd") %>%
  dplyr::mutate(query = gsub("AF-", "", query),
                query = gsub("-F1-model_V4.pdb", "", query))

# make a make of the genbank accession (gb) to the uniprot accession
gb_to_uniprot <- readr::read_tsv(snakemake@input[['uniprot_acc']]) %>%
  dplyr::rename("genbank" = "From", "uniprot" = "Entry", "uniprot_name" = "Entry Name")

fsk <- dplyr::left_join(fsk, gb_to_uniprot, by = c("query" = "uniprot"))
readr::write_tsv(fsk, snakemake@output[['all_fsk']])
