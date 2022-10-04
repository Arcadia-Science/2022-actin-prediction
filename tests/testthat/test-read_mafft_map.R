test_that("check dataframe dimensions are accurate.", {
  df <- read_mafft_map(test_path("O24426_CHLRE.fasta.map"))
  expect_equal(ncol(df), 3)
  expect_equal(nrow(df), 380)
})
