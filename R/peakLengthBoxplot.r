#' Produces box-and-whisker plot of peak lengths within any genomic interval.
#'
#' Makes boxplots of peak lengths to show how lengths of individual peaks are distributed within any genomic interval (e.g., 0-500 bp upstream of nearest gene for all genes throughout the genome).
#'
#' @param organism Object name assigned from readGFF() command.
#' @param start Lower bound of upstream extension.
#' @param end Upper bound of upstream extension.
#'
#' @return Creates boxplots showing how lengths of peaks are distributed within any given genomic interval.  Also, creates character vector composed of individual peak lengths.
#'
#' @examples
#' rat <- readGFF("ftp://ftp.ensembl.org/pub/release-84/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.84.gtf.gz")
#' allpeaks <- system.file("extdata", "totalpeaksfile.txt", package="geneXtendeR")
#' peaksInput(allpeaks)
#' peakLengthBoxplot(rat, 0, 500)
#'
#'
#' @export
peakLengthBoxplot <- function(organism, start, end) {
  a <- distinct(organism, start, end)
  b <- as.numeric(a$V3) - as.numeric(a$V2)
  boxplot(b, ylab = "peak length (bp)")
  b
}