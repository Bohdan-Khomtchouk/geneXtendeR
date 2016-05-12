#' Annotate peaks file.
#'
#' Annotate peaks file with gene information based on optimally chosen geneXtendeR upstream extension file.
#' 
#' @param organism Name of organism (refer to gtf R package).
#' @param extension Desired upstream extension.
#'
#' @return A file of peaks annotated with gene ID, gene name, and gene-to-peak genomic distance (in bp).
#'
#' @examples
#' \dontrun{annotate(rat, 2500)}
#'
#' @export
annotate <- function(organism, extension) {
		onegxFile <- sprintf("geneXtender_gtf_%s.txt", extension)
		options(warn = -1)
		onecmd <- sprintf("bedtools closest -d -a peaks.txt -b %s", onegxFile)
		onecmd.bedtools.closest.output <- system(onecmd, intern = TRUE, ignore.stderr = TRUE)
		write.table(read.table(text = onecmd.bedtools.closest.output), 
		            file = sprintf("peaks_annotated_%s.txt", extension),
		            quote = FALSE, 
		            sep = "\t", 
		            row.names = FALSE,
		            col.names = c("Chromosome", "Peak Start", "Peak End", "Chromosome", "Gene Start", "Gene End", "Gene ID", "Gene Name", "Distance-of-Gene-to-Nearest-Peak")
		            )
	}