test_that("diagonal matrixis formatted correctly and contains correct information about calculating PID", {
  source("../R/pairwise_identity.R")
  pid_mat <- calculate_pid_matrix("confident_alignment.fasta")
  expect_equal(table(is.na(pid_mat))[[1]], 120) # 120 FALSE
  expect_equal(table(!is.na(pid_mat))[[1]], 105) # 105 FALSE
  expect_equal(colnames(pid_mat)[1:(ncol(pid_mat) - 1)], rownames(pid_mat)[2:nrow(pid_mat)]) # staggered off by one because no diagonal
  expect_equal(max(pid_mat, na.rm = T) <= 1, T) # no value should exceed 1
  expect_equal(min(pid_mat, na.rm = T) >= 0, T) # no value should be less than 0

})

test_that("pid_matrix_to_df produces a df of the correct length when there's no filtering", {
  source("../R/pairwise_identity.R")
  pid_mat <- calculate_pid_matrix("confident_alignment.fasta")
  pid_df <- pid_matrix_to_df(pid_mat, filter_to_query_protein = F)
  expect_equal(nrow(pid_df), 120) # after removing NAs, we should have 120 rows left over
})

test_that("pid_matrix_to_df produces a df of the correct length when there is filtering", {
  source("../R/pairwise_identity.R")
  pid_mat <- calculate_pid_matrix("confident_alignment.fasta")
  pid_df <- pid_matrix_to_df(pid_mat, filter_to_query_protein = T)
  expect_equal(nrow(pid_df), 15) # filtering only to the query protein, we should have 15 matches
  expect_equal(length(unique(pid_df$prot2)), 1) # only one protein should be represented in the prot2 column
})

test_that("all pid functions are integrated together correctly", {
  source("../R/pairwise_identity.R")
  msa_file <- "confident_alignment.fasta"
  function_output <- read_msa_and_plot_pid(msa_file, filter_to_query_protein = T, tsv = NULL, pdf = NULL)

  pid_mat <- calculate_pid_matrix(msa_file)
  pid_df <- pid_matrix_to_df(pid_mat, filter_to_query_protein = T, tsv = NULL)
  pid_df_for_plot <- pid_matrix_to_df(pid_mat, filter_to_query_protein = F)
  pid_plt <- plot_pid_df(pid_df = pid_df_for_plot, pdf = NULL)
  nonfunction_output <- pid_plt

  expect_equal(function_output, nonfunction_output)
})
