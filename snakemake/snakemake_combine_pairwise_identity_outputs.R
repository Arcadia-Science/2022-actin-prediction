library(readr)
library(dplyr)
library(purrr)

pid_files <- unlist(snakemake@input[['tsv']])

pid <- pid_files %>%
  purrr::map_dfr(readr::read_tsv, show_col_types = F)

readr::write_tsv(pid, file = snakemake@output[['pid']])

avg_pid <- pid %>%
  dplyr::group_by(prot2) %>%
  dplyr::summarise(avg_fid = mean(fid),
                   avg_pid = mean(pid))

readr::write_tsv(avg_pid, file = snakemake@output[['avg_pid']])
