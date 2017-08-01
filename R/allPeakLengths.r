#' Produces box-and-whisker plot showing distribution of peak lengths across a peaks input file.
#'
#' Makes boxplots of all peak lengths (within a peaks input file) to show how lengths of individual peaks are distributed across the entire peak set.
#'
#' @param filename Name of peaks input file.
#'
#' @return Returns a box-and-whisker plot of peak length distribution across a peaks file.
#'
#' @examples
#' myfile <- system.file("extdata", "somepeaksfile.txt", package="geneXtendeR")
#' allPeakLengths(myfile)
#'
#' @import data.table
#' 
#' @export
allPeakLengths <- function(filename) {
  oopts = options(warn=-1)
  on.exit(options(oopts))
  file.input <- data.table::fread(filename)
  file.input$diff <- file.input$end - file.input$start
  boxplot(file.input$diff, horizontal = TRUE, xlab = "peak length (bp)")
}