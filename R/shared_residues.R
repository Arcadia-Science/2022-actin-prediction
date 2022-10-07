#' combine_alignment_and_feature
#'
#' Combine an alignment and feature file into a single tibble.
#' Requires a tibble read in from read_mafft_map (alignment) and a tibble read in with read_feature_df (feature file).
#' If the file argument is specified, the function will output a tsv of amino acid identities for the query and reference protein at feature residues.
#'
#' @param mafft_map
#' @param feature_df
#' @param file
#'
#' @return A tibble.
#' @export
#'
#' @examples
combine_alignment_and_feature <- function(mafft_map, feature_df, file = NULL){
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop(
      "Package \"dplyr\" must be installed to use this function.",
      call. = FALSE
    )
  }

  df <- feature_df %>%
    dplyr::left_join(mafft_map, by = c("position_reference" = "position_reference_alignment")) %>%
    dplyr::select(feature, position_reference, letter_reference, position_query, letter_query)

  if (!is.null(file)) {
      if (!requireNamespace("readr", quietly = TRUE)) {
        stop(
          "Package \"readr\" must be installed if file is not NULL.",
          call. = FALSE
        )
      }
    readr::write_tsv(df, file = file)
  }

  return(df)
}

#' calculate_shared_residues
#'
#' Calculate the number and fraction of residues with shared identity between a reference protein and a query protein.
#' Requires a tibble generated from combine_alignment_and_feature
#' If the file argument is specified, the function will output a tsv that summarizes the number of compared residues, the number of matching residues, and the fraction of matching residues.
#'
#' @param combined_feature_and_alignment
#' @param file
#'
#' @return A tibble.
#' @export
#'
#' @examples
calculate_shared_residues <- function(combined_feature_and_alignment, file = NULL){
  # guard against missing package installation
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop(
      "Package \"dplyr\" must be installed to use this function.",
      call. = FALSE
    )
  }
  
  # shared residue summary functionality
  df_summary <- combined_feature_and_alignment %>%
    dplyr::mutate(query_matches_reference = ifelse(letter_reference == letter_query, TRUE, FALSE)) %>%
    dplyr::group_by(feature) %>%
    dplyr::summarise(feature_count = dplyr::n(),
                     num_matching = sum(query_matches_reference),
                     fraction_matching = sum(query_matches_reference) / sum(feature_count))
  
  # write file
  if (!is.null(file)) {
      if (!requireNamespace("readr", quietly = TRUE)) {
        stop(
          "Package \"readr\" must be installed if file is not NULL.",
          call. = FALSE
        )
      }
    readr::write_tsv(df_summary, file = file)
  }

  return(df_summary)
}