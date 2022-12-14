test_that("check that join didn't introduce NAs or otherwise cull information", {
  source("../R/shared_residues.R")
  source("../R/read_mafft_map.R")
  source("../R/read_feature_csv.R")
  map <- read_mafft_map("O24426_CHLRE.fasta.map")
  feat <- read_feature_csv("../inputs/protein_features/atp_binding.csv")
  df <- combine_alignment_and_feature(mafft_map = map, feature_df = feat)
  expect_equal(nrow(feat), nrow(df))
  expect_equal(nrow(df[complete.cases(df), ]), nrow(df))
})
