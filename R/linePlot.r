#' Produces line plots.
#'
#' Makes differential line plots showing the differences in the number of genes under peaks at consecutive upstream extension levels.
#'
#' @param organism Object name assigned from readGFF() command.
#' @param start Lower bound of upstream extension.
#' @param end Upper bound of upstream extension.
#' @param by Interval between consecutive extensions.
#'
#' @return Creates differential line plots.
#'
#' @examples
#' rat <- readGFF("ftp://ftp.ensembl.org/pub/release-84/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.84.gtf.gz")
#' fpath <- system.file("extdata", "somepeaksfile.txt", package="geneXtendeR")
#' peaksInput(fpath)
#' linePlot(rat, 1000, 3000, 100)
#'
#'
#' @export
linePlot <- function(organism, start, end, by) {
 if(!file.exists("peaks.txt")){
   message("Please run peaksInput() function first!  See ?peaksInput for more information")
 } else {
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
	par(mar = c(8.1,4.1,2.1,2.1))
	plot(differences, type = "o", col = "blue", xaxt = "n", ylab = "differences", xlab = "")
	axis(1, at = 1:length(differences), labels = xDeltas, las = 3)
	mtext(side = 1, "Genomic region (bp)", line = 6.6)
	}
}
