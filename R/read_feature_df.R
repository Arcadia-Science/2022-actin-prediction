#' read_feature_df
#'
#' Read in a hand-curated feature df.
#' The first column is the numerical position in the amino acid sequence (1-based, integer).
#' The second column is the letter of the amino acid residue.
#' The third column is the feature annotation (any column name).
#' If the feature is observed at a given residue, it's name should occur in the third column.
#' Each csv file should only include annotations for a single feature type (e.g. atp_binding).
#'
#' @param feature_df_file
#'
#' @return A tibble.
#' @export
#'
#' @examples
#' read_feature_df("inputs/atp_binding_doi_10.1073_pnas.1807028115.csv")
read_feature_df <- function(feature_df_file){
  feature_df <- readr::read_csv(feature_df_file,
                                trim_ws = T,           # just in case clean up white spaces so integers and NAs read in properly
                                col_names = c("position_reference", "letter_reference", "feature"), # standardize column names
                                skip = 1,              # skip built in header row
                                col_types = "dcc") %>% # explicitly set column data types
    dplyr::mutate(letter_reference = toupper(letter_reference)) %>% # make sure letters are upper case to match with mafft mapout
    dplyr::filter(!is.na(feature)) # filter to positions annotated with feature
  return(feature_df)
}
