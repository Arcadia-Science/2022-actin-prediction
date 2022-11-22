library(dplyr)
library(readr)
library(purrr)
results <- snakemake@input[['results']] %>%
  purrr::map_dfr(read_tsv)

readr::write_tsv(results, snakemake@output[['tsv']])
