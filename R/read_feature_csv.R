#' Read in hand-curated feature csv file to tibble.
#'
#' Read in a hand-curated feature csv file to a tibble.
#' The first column is the numerical position in the amino acid sequence (1-based, integer).
#' The second column is the letter of the amino acid residue.
#' The third column is the feature annotation (any column name).
#' If the feature is observed at a given residue, it's name should occur in the third column.
#' Each csv file should only include annotations for a single feature type (e.g. atp_binding).
#'
#' @param feature_csv Path to csv file containing feature information.
#'
#' @return A tibble.
#' @export
#'
#' @examples
#' read_feature_csv("inputs/protein_features/atp_binding.csv")
read_feature_csv <- function(feature_csv){
  # functionality to read in feature data frame
  feature_df <- readr::read_csv(feature_csv,
                                trim_ws = T,           # clean up white spaces so integers and NAs read in properly
                                col_names = c("position_reference", "letter_reference", "feature"), # standardize column names
                                skip = 1,              # skip built in header row
                                col_types = "dcc") %>% # explicitly set column data types
    dplyr::mutate(letter_reference = toupper(letter_reference)) %>% # make sure letters are upper case to match with mafft mapout
    dplyr::filter(!is.na(feature)) # filter to positions annotated with feature
  return(feature_df)
}
