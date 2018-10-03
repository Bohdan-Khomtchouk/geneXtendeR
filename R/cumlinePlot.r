#' Produces cumulative line plots.
#'
#' Makes cumulative differential line plots showing the cumulative sums of the number of genes under peaks at consecutive upstream extension levels.
#'
#' @param organism Object name assigned from readGFF() command.
#' @param start Lower bound of upstream extension.
#' @param end Upper bound of upstream extension.
#' @param by Interval between consecutive extensions.
#'
#' @return Creates cumulative differential line plots.
#'
#' @examples
#' library(rtracklayer)
#' rat <- readGFF("ftp://ftp.ensembl.org/pub/release-84/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.84.gtf.gz")
#' fpath <- system.file("extdata", "somepeaksfile.txt", package="geneXtendeR")
#' peaksInput(fpath)
#' cumlinePlot(rat, 1000, 3000, 100)
#'
#'
#' @export
cumlinePlot <- function(organism, start, end, by) {
  if(!file.exists("peaks.txt")){
    message("Please run peaksInput() function first!  See ?peaksInput for more information")
  } else {
    oopts = options(warn=-1)
    on.exit(options(oopts))
    
    sapply(seq(start, end, by), .geneXtender, organism, FALSE)
    xlabs <- as.character(seq(start, end, by))
    xDeltas <- vapply(seq_along(xlabs)[-1], function(i) paste(xlabs[(i-1):i], collapse = "-"), "")
    gxFiles <- sprintf("geneXtender_gtf_%s.bed", seq(start, end, by))
    
    
    run <- function(f1, f2, num) {
      .C("extractnumber", f1, f2, num)[[3]]
    }
    
    num = 0
    numvec <- numeric()
    for (x in gxFiles) {
      numbers <- run(f1 = "peaks.txt", f2 = x, as.integer(num))
      numvec <- append(numvec, numbers)
    }
    
    for (y in numvec) {
      diffs <- sapply(1:length(numvec), function(i) {numvec[i+1] - numvec[i]})
    }
    
    differences <- diffs[!is.na(diffs)]
    cumulative_differences <- cumsum(differences)
    par(mar = c(8.1,4.1,2.1,2.1))
    plot(cumulative_differences, type = "o", col = "red", xaxt = "n", ylab = "cumulative differences", xlab = "")
    axis(1, at = 1:length(cumulative_differences), labels = xDeltas, las = 3)
    mtext(side = 1, "Genomic region (bp)", line = 6.5)
    
  }
}
