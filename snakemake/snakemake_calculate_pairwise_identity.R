source("R/pairwise_identity.R")
library(tidyverse)

read_msa_and_plot_pid(msa_file = snakemake@input[['msa']],
                      filter_to_query_protein = T, 
                      tsv = snakemake@output[['tsv']],
                      pdf = snakemake@output[['pdf']])
