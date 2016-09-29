#' Sample peaks list to be used as input to geneXtendeR
#'
#' A dataset containing the chromosome number, start and stop positions
#' of ChIP-seq peaks along the Rattus norvegicus genome (rn6 assembly).
#' A dataset like this may be used as input to the peaksInput() command, 
#' which will sort the dataset by chromosome number and start position.
#'
#' @docType data
#'
#' @usage data(samplepeaksinput)
#' 
#' @return Demonstrates a sample peaks file used as input.
#'
#' @examples 
#' head(samplepeaksinput)
#' tail(samplepeaksinput)
#'
#' @format A data frame with 25089 rows and 3 variables:
#' \describe{
#' \item{chr}{Chromosome number}
#' \item{start}{Peak start position [in units of base pairs]}
#' \item{end}{Peak end position [in units of base pairs]}
#' }
#'
#' @keywords datasets
"samplepeaksinput"
