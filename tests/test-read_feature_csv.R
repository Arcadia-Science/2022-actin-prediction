test_that("check that no NAs sneak in to the feature df", {
  source("../R/read_feature_df.R")
  df <- read_feature_csv("../inputs/protein_features/atp_binding.csv")
  expect_equal(nrow(df[complete.cases(df), ]), nrow(df))
})
