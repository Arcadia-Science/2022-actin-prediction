library(dplyr)
library(readr)
library(purrr)
results <- snakemake@input[['results']] %>%
  map_dfr(read_tsv)

write_tsv(results, snakemake@output[['tsv']])
