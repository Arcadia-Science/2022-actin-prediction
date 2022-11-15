library(readr)

results <- read_tsv("outputs/foldseek/uniprot_accessions/results.tsv")

# check if the output dir exists, and if not, create it:
if(dir.exists(snakemake@output[['outdir']]) == FALSE){
  dir.create(snakemake@output[['outdir']])
}

for(entry in results$Entry){
  file.create(paste0(snakemake@output[['outdir']], "/", entry, ".txt"))
}
