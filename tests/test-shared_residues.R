test_that("check that giardia residue calculations match expectations", {
  source("../R/shared_residues.R")
  df <- readr::read_tsv("P51775_ACT_GIAIN-longitudinal_actin_contact_full.tsv")
  df_summary <- calculate_shared_residues(df)
  test_summary <- readr::read_tsv("P51775_ACT_GIAIN-longitudinal_actin_contact_summary.tsv")
  expect_equal(df_summary$num_matching, test_summary$num_matching)
})

test_that("check that bovine residue calculations match expectations", {
  source("../R/shared_residues.R")
  df <- readr::read_tsv("P63258_ACTG_BOVIN-atp_binding_full.tsv")
  df_summary <- calculate_shared_residues(df)
  test_summary <- readr::read_tsv("P63258_ACTG_BOVIN-atp_binding_summary.tsv")
  expect_equal(df_summary$num_matching, test_summary$num_matching)
})
