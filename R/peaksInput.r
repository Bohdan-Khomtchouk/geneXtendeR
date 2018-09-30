#' Preprocesses a peaks input file.
#'
#' Takes your tab-delimited 3-column (chromosome number, peak start, and peak end) input file (see ?samplepeaksinput) consisting of peaks called from a peak caller (e.g., MACS2 or SICER) and sorts the file by chromosome and start position, thereby creating a preprocessed file for downstream geneXtendeR analysis.  This file (called "peaks.txt") is a preprocessed file of the original input and is deposited in the user's working directory and used for the remainder of the analysis.  In this "peaks.txt" file, peaks are sorted by chromosome number and start position, where the X chromosome is designated by the integer 100, the Y chromosome by the integer 200, and the mitochondrial chromosome by the integer 300.
#'
#' @param filename Name of file containing peaks that have been generated from a peak caller (e.g., MACS2, SICER).  See ?samplepeaksinput for an example of such an input file.
#'
#' @return Returns a formatted file (called "peaks.txt") that has been preprocessed in preparation for usage with barChart(), linePlot(), distinct(), and other downstream commands and deposited in the user's working directory.
#'
#' @examples
#' ?samplepeaksinput  #Documentation of some exemplary sample input
#' data(samplepeaksinput)
#' head(samplepeaksinput)
#' tail(samplepeaksinput)
#' fpath <- system.file("extdata", "somepeaksfile.txt", package="geneXtendeR")
#' peaksInput(fpath)
#'
#' @importFrom data.table fread
#' @import rtracklayer
#' 
#' @export
peaksInput <- function(filename) {
  file.input <- fread(filename)
  if (ncol(file.input) != 3) {
    stop("Must be a 3-column tab-demimited file. (see ?samplepeaksinput)")
  }
	file.input[chr == "chr1", chr := '1']
	file.input[chr == "chr2", chr := '2']
	file.input[chr == "chr3", chr := '3']
	file.input[chr == "chr4", chr := '4']
	file.input[chr == "chr5", chr := '5']
	file.input[chr == "chr6", chr := '6']
	file.input[chr == "chr7", chr := '7']
	file.input[chr == "chr8", chr := '8']
	file.input[chr == "chr9", chr := '9']
	file.input[chr == "chr10", chr := '10']
	file.input[chr == "chr11", chr := '11']
	file.input[chr == "chr12", chr := '12']
	file.input[chr == "chr13", chr := '13']
	file.input[chr == "chr14", chr := '14']
	file.input[chr == "chr15", chr := '15']
	file.input[chr == "chr16", chr := '16']
	file.input[chr == "chr17", chr := '17']
	file.input[chr == "chr18", chr := '18']
	file.input[chr == "chr19", chr := '19']
	file.input[chr == "chr20", chr := '20']
	file.input[chr == "chr21", chr := '21']
	file.input[chr == "chr22", chr := '22']
	file.input[chr == "chr23", chr := '23']
	file.input[chr == "chr24", chr := '24']
	file.input[chr == "chr25", chr := '25']
	file.input[chr == "chr26", chr := '26']
	file.input[chr == "chr27", chr := '27']
	file.input[chr == "chr28", chr := '28']
	file.input[chr == "chr29", chr := '29']
	file.input[chr == "chr30", chr := '30']
	file.input[chr == "chr31", chr := '31']
	file.input[chr == "chr32", chr := '32']
	file.input[chr == "chr33", chr := '33']
	file.input[chr == "chr34", chr := '34']
	file.input[chr == "chr35", chr := '35']
	file.input[chr == "chr36", chr := '36']
	file.input[chr == "chr37", chr := '37']
	file.input[chr == "chr38", chr := '38']
	file.input[chr == "chr39", chr := '39']
	file.input[chr == "chr40", chr := '40']
	file.input[chr == "chr41", chr := '41']
	file.input[chr == "chr42", chr := '42']
	file.input[chr == "chr43", chr := '43']
  file.input[chr == "X", chr := '100']
  file.input[chr == "chrX", chr := '100']
  file.input[chr == "Y", chr := '200']
  file.input[chr == "chrY", chr := '100']
  file.input[chr == "MT", chr := '300']
  file.input[chr == "chrMT", chr := '300']
  file.input[chr == "MtDNA", chr := '300']
  file.input[chr == "chrMtDNA", chr := '300']
  file.input[chr == "M", chr := '300']
  file.input[chr == "chrM", chr := '300']
  file.input[chr == "Mito", chr := '300']
  file.input[chr == "chrMito", chr := '300']
  file.sorted <- dplyr::arrange(file.input, as.numeric(chr), as.numeric(start))
  write.table(na.omit(file.sorted), file = "peaks.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}
