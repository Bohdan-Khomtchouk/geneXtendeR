#' Peaks list ordered by chromosome number
#'
#' A dataset containing the start and stop positions
#' of ChIP-seq peaks along a genome
#'
#' @docType data
#'
#' @usage data(sample_peaks_input)
#' 
#' @examples 
#' head(sample_peaks_input)
#' tail(sample_peaks_input)
#'
#' @format A data frame with 9799 rows and 3 variables:
#' \describe{
#' \item{V1}{Chromosome number}
#' \item{V2}{Peak start position [in units of base pairs]}
#' \item{V3}{Peak end position [in units of base pairs]}
#' }
#'
#' @keywords datasets
"sample_peaks_input"
