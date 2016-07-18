#' Peaks list ordered by chromosome number and start position
#'
#' A dataset containing the chromosome number, start and stop positions
#' of ChIP-seq peaks along a genome
#'
#' @docType data
#'
#' @usage data(samplepeaksinput)
#' 
#' @examples 
#' head(samplepeaksinput)
#' tail(samplepeaksinput)
#'
#' @format A data frame with 25088 rows and 3 variables:
#' \describe{
#' \item{V1}{Chromosome number}
#' \item{V2}{Peak start position [in units of base pairs]}
#' \item{V3}{Peak end position [in units of base pairs]}
#' }
#'
#' @keywords datasets
"samplepeaksinput"
