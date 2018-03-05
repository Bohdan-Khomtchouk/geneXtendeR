#' Annotate peak file based on gene.
#'
#' Annotate a user's peaks file (which has been preprocessed with the peaksInput() command) with gene information based on optimally chosen geneXtendeR upstream extension file and compresses the annotations based on genes.  This command requires a preprocessed "peaks.txt" file (generated using peaksInput()) to be present in the user's working directory, otherwise the user is prompted to rerun the peaksInput() command in order to regenerate it.
#'
#' @param organism Object name assigned from readGFF() command.
#' @param extension Desired upstream extension.
#' 
#' @return The gene coordinates are extended by `extension` at the 5-prime end, and by 500 bp at the 3-prime end.  The peaks file is then overlayed on these new gene coordinates, producing a file of peaks annotated with gene ID, gene name, gene location, mean and standard deviation of peaks-to-genes, number of peaks-to-genes, and peak-to-gene genomic distance (in bp).  Distance is calculated between 5-prime end of gene and 3-prime end of peak.
#' @return A data.table formatted version of the gene-annotated file for checking or further calculations.
#' 
#' @return (From annotate.r) The gene coordinates are extended by `extension` at the 5-prime end, and by 500 bp at the 3-prime end.  The peaks file is then overlayed on these new gene coordinates, producing a file of peaks annotated with gene ID, gene name, and gene-to-peak genomic distance (in bp).  Distance is calculated between 5-prime end of gene and 3-prime end of peak.
#'
#' @examples
#' rat <- readGFF("ftp://ftp.ensembl.org/pub/release-84/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.84.gtf.gz")
#' fpath <- system.file("extdata", "somepeaksfile.txt", package="geneXtendeR")
#' peaksInput(fpath)
#' gene_annotate(rat, 3400)
#' 
#' @useDynLib geneXtendeR, .registration = TRUE
#' @import data.table
#'
#' @export
gene_annotate <- function(organism, extension) {
  PeaktoGeneAnno <- annotate(organism, extension)
  
  GtPA <- PeaktoGeneAnno[,
                         .('Peaks-on-Gene-Body' = sum(`Distance-of-Gene-to-Nearest-Peak` == 0),
                           mean = mean(`Distance-of-Gene-to-Nearest-Peak`),
                           sd = sd(`Distance-of-Gene-to-Nearest-Peak`),
                           .N),
                         by=.(`Chromosome`, `Gene-Start`, `Gene-End`, `Gene-ID`, `Gene-Name`)]
  data.table::setnames(GtPA, c("N", "mean"), c("Number-of-Peaks-Associated-with-Gene", "Mean-Distance-of-Gene-to-Nearest-Peaks"))
  
  GtPA <- dplyr::arrange(GtPA,desc(`Peaks-on-Gene-Body`), desc(`Number-of-Peaks-Associated-with-Gene`))
  
  data.table::fwrite(
    GtPA,
    file = sprintf("gene_annotated_%s.txt", extension),
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )
  return(GtPA[])
}