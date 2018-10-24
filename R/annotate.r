#' Annotate peaks file.
#'
#' Annotate a user's peaks file (which has been preprocessed with the peaksInput() command) with gene information based on optimally chosen geneXtendeR upstream extension file.  This command requires a preprocessed "peaks.txt" file (generated using peaksInput()) to be present in the user's working directory, otherwise the user is prompted to rerun the peaksInput() command in order to regenerate it.
#' 
#' @param organism Object name assigned from readGFF() command.
#' @param extension Desired upstream extension.
#'
#' @return The gene coordinates are extended by `extension` at the 5-prime end, and by 500 bp at the 3-prime end.  The peaks file is then overlayed on these new gene coordinates, producing a file of peaks annotated with gene ID, gene name, and gene-to-peak genomic distance (in bp).  Distance is calculated between 5-prime end of gene and 3-prime end of peak.
#' @return A data.table formatted version of the annotated file for checking or further calculations.
#'
#' @examples
#' library(rtracklayer)
#' rat <- readGFF("ftp://ftp.ensembl.org/pub/release-84/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.84.gtf.gz")
#' fpath <- system.file("extdata", "somepeaksfile.txt", package="geneXtendeR")
#' peaksInput(fpath)
#' annotate(rat, 2500)
#'
#' @import data.table
#'
#' @useDynLib geneXtendeR, .registration = TRUE
#'
#' @export
annotate <- function(organism, extension) {
  if(!file.exists("peaks.txt")){
    stop("Please run peaksInput() function first!  See ?peaksInput for more information")
  } else {
    oopts = options(warn=-1)
    on.exit(options(oopts))
    
    run3 <- function(f1, f2, peakslist) {
      .C("annotate", f1, f2, peakslist)[[3]]
    }
    
    linelen <- "                                                                                                                                  "
    n <- 500000
    peaksArray<-rep(linelen,n)
    sapply(extension, .geneXtender, organism, FALSE)
    onegxFile <- sprintf("geneXtender_gtf_%s.bed", extension)
    onecmd2 <- run3(f1 = "peaks.txt", f2 = onegxFile, peakslist = peaksArray) 
    onecmd3 <- onecmd2[onecmd2 != linelen]
    
    write.table(
      onecmd3,
      file = sprintf("peaks_annotated_%s.txt", extension),
      sep = "\t",
      row.names = FALSE,
      col.names = paste("Chromosome\t", "Peak-Start\t", "Peak-End\t", "Chromosome\t", "Gene-Start\t", "Gene-End\t", "Gene-ID\t", "Gene-Name\t", "Distance-of-Gene-to-Nearest-Peak"),
      quote = FALSE
    )
    
    rdt <- data.table::transpose(data.table::as.data.table(strsplit(onecmd3, "\t"))) #Formatting
    colnames(rdt) <- c("Chromosome", "Peak-Start", "Peak-End", "rm", "Gene-Start", "Gene-End", "Gene-ID", "Gene-Name", "Distance-of-Gene-to-Nearest-Peak")
    rdt <- rdt[,rm:= NULL]
    for (j in c("Chromosome", "Peak-Start", "Peak-End", "Gene-Start", "Gene-End", "Distance-of-Gene-to-Nearest-Peak")) data.table::set(rdt, j=j, value = as.numeric(rdt[[j]]))
    return(rdt[])
 }      
}