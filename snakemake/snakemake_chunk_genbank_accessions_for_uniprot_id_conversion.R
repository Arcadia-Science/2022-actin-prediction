library(dplyr)
library(readr)

# gb_accs stands for genbank accessions.
# it's solved from the input fasta files, which are named by their genbank protein accessions.
gb_accs <- basename(unlist(snakemake@input[['fastas']]))
gb_accs <- gsub(".fasta", "", gb_accs)

# ngrps is the number of groups to chunk the file into.
# we used 15000 as the chunk size.
# the code following ngrps creates a vector that will be used by dplyr::group_by to write each chunk of 15k to it's own file.
# it deals with potentially uneven sizes from the chunking.

ngrps <- ceiling(length(gb_accs) / 15000)
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
  dplyr::mutate(grp = grp) %>% 
  dplyr::group_by(grp) %>%
  dplyr::group_walk(~ readr::write_tsv(.x, paste0(snakemake@output[['outdir']], "/genbank_accessions_chunk", .y$grp, ".txt"),
                         col_names = F))
