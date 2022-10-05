test_that("check dataframe dimensions are accurate.", {
  source("../R/read_mafft_map.R")
  df <- read_mafft_map("O24426_CHLRE.fasta.map")
  expect_equal(ncol(df), 3)
  expect_equal(nrow(df), 380)
})
