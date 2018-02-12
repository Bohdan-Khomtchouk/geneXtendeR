#' Produces bar charts.
#'
#' Makes bar graphs showing the number of genes under peaks at various upstream extension levels.
#'
#' @param organism Object name assigned from readGFF() command.
#' @param start Lower bound of upstream extension.
#' @param end Upper bound of upstream extension.
#' @param by Interval between consecutive extensions.
#'
#' @return Creates bar charts.
#'
#' @examples
#' rat <- readGFF("ftp://ftp.ensembl.org/pub/release-84/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.84.gtf.gz")
#' fpath <- system.file("extdata", "somepeaksfile.txt", package="geneXtendeR")
#' peaksInput(fpath)
#' barChart(rat, 1000, 3000, 100)
#'
#' @importFrom dplyr filter select arrange
#' @import utils
#' @import graphics
#'
#' @useDynLib geneXtendeR, .registration = TRUE
#'
#'
#' @export
barChart <- function(organism, start, end, by) {
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
	
	barplot(numvec, names.arg = xlabs)
	
	}
}
