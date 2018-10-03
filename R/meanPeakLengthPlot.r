#' Produces line plots of mean (average) peak length within any genomic interval.
#'
#' Makes line plots of mean peak lengths to show the average length of individual peaks within any genomic interval (e.g., 0-500 bp upstream of nearest gene for all genes throughout the genome).
#' 
#' @param organism Object name assigned from readGFF() command.
#' @param start Lower bound of upstream extension.
#' @param end Upper bound of upstream extension.
#' @param by Interval between consecutive extensions.
#'
#' @return Creates mean peak length line plots.
#'
#' @examples
#' library(rtracklayer)
#' rat <- readGFF("ftp://ftp.ensembl.org/pub/release-84/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.84.gtf.gz")
#' allpeaks <- system.file("extdata", "totalpeaksfile.txt", package="geneXtendeR")
#' peaksInput(allpeaks)
#' meanPeakLengthPlot(rat, 0, 10000, 500)
#'
#'
#' @export
meanPeakLengthPlot <- function(organism, start, end, by) {
  x <- seq(start, end, by)
  peak_lengths <- sapply(1:(length(x)-1), function(i) {meanPeakLength(organism, x[i], x[i+1])})
  par(mar = c(7.1,4.1,2.1,2.1))
  plot(peak_lengths, type = "o", col = "blue", xaxt = "n", ylab = "mean peak length (bp)", xlab = "")
  xlabs <- as.character(seq(start, end, by))
  xDeltas <- vapply(seq_along(xlabs)[-1], function(i) paste(xlabs[(i-1):i], collapse = "-"), "")
  axis(1, at = 1:length(peak_lengths), labels = xDeltas, las = 3)
  mtext(side = 1, "Genomic region (bp)", line = 6)
}