#' Calculates mean (average) peak length for any genomic region.
#'
#' Determines the average peak length of all peaks found within some genomic interval (e.g., 0-500 bp upstream of nearest gene for all genes throughout the genome).
#' 
#' 
#' @param organism Object name assigned from readGFF() command.
#' @param start Lower bound of upstream extension.
#' @param end Upper bound of upstream extension.
#'
#' @return A vector composed of a single number representing the average peak length found within a genomic interval.
#'
#' @examples
#' rat <- readGFF("ftp://ftp.ensembl.org/pub/release-84/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.84.gtf.gz")
#' sigpeaks <- system.file("extdata", "significantpeaksfile.txt", package="geneXtendeR")
#' peaksInput(sigpeaks)
#' meanPeakLength(rat, 0, 500)
#'
#' @useDynLib geneXtendeR, .registration = TRUE
#'
#' @export
meanPeakLength <- function(organism, start, end) {
  if(!file.exists("peaks.txt")){
    stop("Please run peaksInput() function first!  See ?peaksInput for more information")
  } else {
    
    run2 <- function(f1, f2, peakslist) {
      .C("extractpeaks", f1, f2, peakslist)[[3]]
    }
    
    sapply(c(start, end), .geneXtender, rat, FALSE)
    twogxFiles <- sprintf("geneXtender_gtf_%s.bed", c(start, end))
    linelen <- ""  
    n <- 500000
    peaksArray <- rep(linelen, n)
    peaksArray2 <- rep(linelen, n)
    cmdtmp1 <- run2(f1 = "peaks.txt", f2 = twogxFiles[[1]], as.character(peaksArray))
    cmdtmp2 <- run2(f1 = "peaks.txt", f2 = twogxFiles[[2]], as.character(peaksArray2))
    cmd1 <- cmdtmp1[cmdtmp1 != linelen]
    cmd2 <- cmdtmp2[cmdtmp2 != linelen]
    m = regexec("^(?:[^\t]+\t){3}", cmd1)
    first3.cmd1 = unlist(regmatches(cmd1, m))
    m = regexec("^(?:[^\t]+\t){3}", cmd2)
    first3.cmd2 = unlist(regmatches(cmd2, m))
    finalList = cmd2[!(first3.cmd2 %in% first3.cmd1)]
    DT <- data.table::as.data.table(do.call("rbind", strsplit(finalList, split = "\t")))
    if (nrow(DT) == 0) return(0)
    return(mean(as.numeric(DT$V3) - as.numeric(DT$V2)))
  }      
}