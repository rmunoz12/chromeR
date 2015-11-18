#' Human chromesome lengths in basepairs
#'
#' A vector containing the lengths of the 22 human autosomes in basepairs.
#'
#' @source \url{https://en.wikipedia.org/wiki/Human_genome}
"chrome_bp"


#' Segments shared IBD (simulated).
#'
#' A dataset containing simulated output from ersa
#' (\url{https://github.com/rmunoz12/ersa}), where each row describes an IBD
#' segment. The variables are as follows:
#'
#' \itemize{
#'  \item id. identifier for each row.
#'  \item result_id. indentifier for each pair of individuals.
#'  \item indv1, indv2. names of individuals.
#'  \item chromosome. chromosomes where each IBD segment is shared.
#'  \item bp_start. starting basepair for each IBD segment.
#'  \item bp_end. ending basepair for each IBD segment.
#'  \item length. length in cM -- values are completely meaningless in this
#'      example data.
#' }
"test_ibd_segments"
