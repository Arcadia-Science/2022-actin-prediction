#' read_mafft_map
#'
#' Read in a file output by mafft when the --mapout flag is used.
#'
#' @param mafft_map_file
#'
#' @return A tibble.
#' @export
#'
#' @examples
#' read_mafft_map("tests/O24426_CHLRE.fasta.map")
read_mafft_map <- function(mafft_map_file){
  # guard against missing package installations
  if (!requireNamespace("readr", quietly = TRUE)) {
    stop(
      "Package \"readr\" must be installed to use this function.",
      call. = FALSE
    )
  }
  
  # read_mafft_map functionality
  mafft_map <- readr::read_csv(mafft_map_file,
                               skip = 2,          # skip mafft-generated metadata rows
                               col_names = c("letter_query", "position_query", "position_reference_alignment"), # reset col names
                               trim_ws = TRUE,    # trim white spaces in csv
                               na = "-",          # set non-aligned from - to NA,
                               col_types = "cdd") # set column types as character, digit, digit
  return(mafft_map)
}
