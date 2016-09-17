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
#' fpath <- system.file("extdata", "somepeaksfile.txt", package="geneXtendeR")
#' peaksInput(fpath)
#'
#' @import data.table
#' @import rtracklayer
#' 
#' @export
peaksInput <- function(filename) {
	options(warn = -1)
	file.input <- fread(filename)
	file.input[V1 == "chr1", V1 := 1]
	file.input[V1 == "chr2", V1 := 2]
	file.input[V1 == "chr3", V1 := 3]
	file.input[V1 == "chr4", V1 := 4]
	file.input[V1 == "chr5", V1 := 5]
	file.input[V1 == "chr6", V1 := 6]
	file.input[V1 == "chr7", V1 := 7]
	file.input[V1 == "chr8", V1 := 8]
	file.input[V1 == "chr9", V1 := 9]
	file.input[V1 == "chr10", V1 := 10]
	file.input[V1 == "chr11", V1 := 11]
	file.input[V1 == "chr12", V1 := 12]
	file.input[V1 == "chr13", V1 := 13]
	file.input[V1 == "chr14", V1 := 14]
	file.input[V1 == "chr15", V1 := 15]
	file.input[V1 == "chr16", V1 := 16]
	file.input[V1 == "chr17", V1 := 17]
	file.input[V1 == "chr18", V1 := 18]
	file.input[V1 == "chr19", V1 := 19]
	file.input[V1 == "chr20", V1 := 20]
	file.input[V1 == "chr21", V1 := 21]
	file.input[V1 == "chr22", V1 := 22]
	file.input[V1 == "chr23", V1 := 23]
	file.input[V1 == "chr24", V1 := 24]
	file.input[V1 == "chr25", V1 := 25]
	file.input[V1 == "chr26", V1 := 26]
	file.input[V1 == "chr27", V1 := 27]
	file.input[V1 == "chr28", V1 := 28]
	file.input[V1 == "chr29", V1 := 29]
	file.input[V1 == "chr30", V1 := 30]
	file.input[V1 == "chr31", V1 := 31]
	file.input[V1 == "chr32", V1 := 32]
	file.input[V1 == "chr33", V1 := 33]
	file.input[V1 == "chr34", V1 := 34]
	file.input[V1 == "chr35", V1 := 35]
	file.input[V1 == "chr36", V1 := 36]
	file.input[V1 == "chr37", V1 := 37]
	file.input[V1 == "chr38", V1 := 38]
	file.input[V1 == "chr39", V1 := 39]
	file.input[V1 == "chr40", V1 := 40]
	file.input[V1 == "chr41", V1 := 41]
	file.input[V1 == "chr42", V1 := 42]
	file.input[V1 == "chr43", V1 := 43]
    file.input[V1 == "X", V1 := 100]
    file.input[V1 == "chrX", V1 := 100]
    file.input[V1 == "Y", V1 := 200]
    file.input[V1 == "chrY", V1 := 100]
    file.input[V1 == "MT", V1 := 300]
    file.input[V1 == "chrMT", V1 := 300]
    file.input[V1 == "MtDNA", V1 := 300]
    file.input[V1 == "chrMtDNA", V1 := 300]
    file.input[V1 == "M", V1 := 300]
    file.input[V1 == "chrM", V1 := 300]
    file.input[V1 == "Mito", V1 := 300]
    file.input[V1 == "chrMito", V1 := 300]
    file.sorted <- dplyr::arrange(file.input, as.numeric(V1), V2)
    write.table(file.sorted, file = "peaks.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
    }