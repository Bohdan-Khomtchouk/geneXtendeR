#' Graphs hotspots of statistically significant peak activity.
#'
#' Makes line plots showing the ratio of statistically significant peaks to the total number of peaks at each genomic interval (e.g., 0-500 bp upstream of every gene in the genome, 500-1000 bp upstream of every gene in the genome, etc.).
#'
#'
#' @param totalpeaksfile Filename in user's working directory (or full path to filename) containing all peaks.
#' @param significantpeaksfile Filename in user's working directory (or full path to filename) containing only significant peaks.
#' @param organism Object name assigned from readGFF() command.
#' @param start Lower bound of upstream extension.
#' @param end Upper bound of upstream extension.
#' @param by Interval between consecutive extensions.
#'
#' @return Line plot showing the ratio of significant to total peaks at each interval across the genome.
#' 
#' @examples
#' rat <- readGFF("ftp://ftp.ensembl.org/pub/release-84/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.84.gtf.gz")
#' allpeaks <- system.file("extdata", "totalpeaksfile.txt", package="geneXtendeR")
#' sigpeaks <- system.file("extdata", "significantpeaksfile.txt", package="geneXtendeR")
#' hotspotPlot(allpeaks, sigpeaks, rat, 0, 10000, 500)
#'
#'
#' @export
hotspotPlot <- function(totalpeaksfile, significantpeaksfile, organism, start, end, by) {
    peaksInput(totalpeaksfile)
  
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
    
    totaldifferences <- diffs[!is.na(diffs)]
    
    #----Now let's repeat for significant peaks----#
    
    peaksInput(significantpeaksfile)
    
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
    
    significantdifferences <- diffs[!is.na(diffs)]
    
    
    psiMetric <- significantdifferences/totaldifferences
    par(mar = c(7.1,4.1,2.1,2.1))
    plot(psiMetric, type = "o", col = "blue", xaxt = "n", ylab = expression(paste("Significant Peaks /Total Peaks ( ", psi, " )")), xlab = "")
    axis(1, at = 1:length(psiMetric), labels = xDeltas, las = 3)
    mtext(side = 1, "Genomic region (bp)", line = 6)
    
  
}
