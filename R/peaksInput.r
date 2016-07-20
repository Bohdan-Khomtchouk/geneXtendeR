#' Formats a peaks input file.
#'
#' Takes your tab-delimited 3-column input file (see Arguments section below) consisting of peaks called from a peak caller (e.g., MACS2 or SICER) and sorts the file by chromosome and start position to create a formatted file for downstream geneXtendeR analysis.  For sorting purposes, the X chromosome is designated by the integer 100, the Y chromosome by the integer 200, and the mitochondrial chromosome by the integer 300.
#'
#' @param filename Name of file containing peaks.  This tab-delimited input file may only contain 3 columns: chromosome number, peak start, and peak end.  See ?samplepeaksinput for more details.
#'
#' @return Returns a formatted file to be used with barChart(), linePlot(), and other downstream commands.
#'
#' @examples
#' ?samplepeaksinput  #Documentation of some exemplary sample input
#' data(samplepeaksinput)
#' head(samplepeaksinput)
#' tail(samplepeaksinput)
#' \dontrun{peaksInput("somepeaksfile.txt")}
#'
#' @import data.table
#' 
#' @export
peaksInput <- function(filename) {
	options(warn = -1)
	file.input <- fread(filename)
    file.input[V1 == "X", V1 := 100]
    file.input[V1 == "Y", V1 := 200]
    file.input[V1 == "MT", V1 := 300]
    file.input[V1 == "MtDNA", V1 := 300]
    file.input[V1 == "M", V1 := 300]
    file.input[V1 == "Mito", V1 := 300]
    file.sorted <- dplyr::arrange(file.input, as.numeric(V1), V2)
    write.table(file.sorted, file = "peaks.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
    }