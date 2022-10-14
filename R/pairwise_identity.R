#' Calculate pairwise identity matrix for a FASTA multiple sequence alignment.
#'
#' @param msa_file Path to a multiple sequence alignment file in fasta format.
#'
#' @return
#' @export
#'
#' @examples
calculate_pid_matrix <- function(msa_file){
  msa <- bio3d::read.fasta(file = msa_file) # read in the MSA
  mat <- bio3d::seqidentity(msa)            # calculate pairwise identity
  mat[lower.tri(mat, diag = TRUE)] <- NA    # fill the lower triangle with NAs 
  mat <- mat[-nrow(mat), -1]                # remove the diagonal
  return(mat)
}

#' Reformat a pairwise identity matrix into a long data.frame.
#'
#' @param pid_mat A matrix produced by calculate_pid_matrix.
#' @param filter_to_query_protein Boolean. Filter to the query protein sequence or return all pairwise identity values.
#' @param tsv Path to output a tsv file of long format pairwise identity values.
#'
#' @return A data.frame.
#' @export
#'
#' @examples
pid_matrix_to_df <- function(pid_mat, filter_to_query_protein = T, tsv = NULL){
  # reformat the pid matrix from wide format to long format 
  df <- pid_mat %>%
    as.data.frame() %>%
    tibble::rownames_to_column("seq1") %>%
    tidyr::pivot_longer(cols = -seq1, names_to = "seq2", values_to = "fid") %>%
    dplyr::filter(!is.na(fid)) %>%
    dplyr::mutate(seq1 = forcats::fct_inorder(seq1),
                  seq2 = forcats::fct_inorder(seq2),
                  pid = fid * 100)
  
  # if a character string representing a protein name is provided, filter to this protein value
  if(filter_to_query_protein == T){
    # This is somewhat brittle but works in the context of our workflow.
    # This approach automatically selects which protein to filter on based on the column name of the last column of the pid_mat.
    # Because we concatenate our confident actin sequences with our query protein, our query protein will always be the last one in the multiple sequence alignment and therefore represented by the last column in the pid_mat.
    # This is better than asking to supply to protein name:
    # In the context of a workflow, a user might encode their query_protein filename, therefore the wildcard, differently than the fasta sequence.
    # This would break the implementation if we filtered on protein name using the wildcard value.
    
    query_protein <- colnames(pid_mat)[ncol(pid_mat)] # extract column name of last column of pid_mat...
    df <- df %>% 
      dplyr::filter(seq2 == query_protein)            # ... and use it to filter the df
  }
  
  # if a file path is provided for the tsv argument, write the data frame to a tsv file
  if(!is.null(tsv)){
    readr::write_tsv(df, file = tsv)
  }
  
  return(df)
} 

#' Use ggplot2 to plot pairwise identity information.
#'
#' @param pid_df Data.frame produced by pid_matrix_to_df with filter_to_query_protein = F.
#' @param pdf Path to output a pdf file of the plot.
#'
#' @return A ggplot2 object.
#' @export
#'
#' @examples
plot_pid_df <- function(pid_df, pdf = NULL){
  plt <- ggplot2::ggplot(pid_df, ggplot2::aes(x = seq1, y = seq2, fill = pid, label = round(pid))) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(color = "white", size = 3) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.title = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 90, hjust = -0.01)) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::scale_fill_viridis_c() +
    ggplot2::coord_equal()
  
  # if a file path is given for pdf, save the image to a pdf
  if(!is.null(pdf)){
    ggplot2::ggsave(filename = pdf, plot = plt)
  } else {
    return(plt)
  }
}

#' Read in a MSA fasta file, calculate pairwise identity, and plot.
#'
#' @param msa_file Path to a multiple sequence alignment file in fasta format.
#' @param filter_to_query_protein  Boolean. Filter to the query protein sequence or return all pairwise identity values.
#' @param tsv Path to output a tsv file of long format pairwise identity values.
#' @param pdf Path to output a pdf file of the plot.
#'
#' @return A ggplot2 object.
#' @export
#'
#' @examples
read_msa_and_plot_pid <- function(msa_file, filter_to_query_protein = T, tsv = NULL, pdf = NULL){
  pid_mat <- calculate_pid_matrix(msa_file)
  pid_df <- pid_matrix_to_df(pid_mat, filter_to_query_protein = filter_to_query_protein, tsv = tsv)
  pid_df_for_plot <- pid_matrix_to_df(pid_mat, filter_to_query_protein = F)
  pid_plt <- plot_pid_df(pid_df = pid_df_for_plot, pdf = pdf)
  return(pid_plt)
}