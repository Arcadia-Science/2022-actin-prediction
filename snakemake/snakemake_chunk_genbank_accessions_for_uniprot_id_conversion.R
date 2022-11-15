library(dplyr)
library(readr)

gb_accs <- basename(unlist(snakemake@input[['fastas']]))
gb_accs <- gsub(".fasta", "", gb_accs)

print(gb_accs)
ngrps <- ceiling(length(gb_accs) / 10000)
grp <- rep(1:ngrps, times = length(gb_accs)/ngrps)

if(length(grp) != length(gb_accs)) {
  difference <-  length(gb_accs) - length(grp)
  if(difference > 0){
    # if the difference is positive, add a vector of 1 to the end of grp
    add <- rep(1, times = difference)
    grp <- c(grp, add)
  } else {
    # if the difference is negative, remove until equal
    grp <- grp[1:(length(grp) + difference)]
  }
}

# stop if the length of the grp vector isn't the same as the number of genbank accessions
all.equal(length(grp), length(gb_accs))

# check if the output dir exists, and if not, create it:
if(dir.exists(snakemake@output[['outdir']]) == FALSE){
  dir.create(snakemake@output[['outdir']])
}
# group by the grp variable and write to output
data.frame(acc = gb_accs) %>% 
  mutate(grp = grp) %>% 
  group_by(grp) %>%
  group_walk(~ write_tsv(.x, paste0(snakemake@output[['outdir']], "/genbank_accessions_chunk", .y$grp, ".txt"),
                         col_names = F))
