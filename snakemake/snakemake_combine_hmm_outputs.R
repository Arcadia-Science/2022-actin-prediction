library(readr)
library(dplyr)
library(purrr)
source(snakemake@input[['rhmmer']])

hmm_files <- unlist(snakemake@input[['tbl']])

hmm <- hmm_files %>%
  purrr::map_dfr(rhmmer::read_tblout)

readr::write_tsv(hmm, snakemake@output[['all_hmm']])
