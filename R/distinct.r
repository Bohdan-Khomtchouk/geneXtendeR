#' Finds unique genes under peaks.
#'
#' Determines what genes directly under peaks are actually unique between two different upstream extension levels.
#' 
#' V1-V3 denote the chromosome/start/end positions of the peaks, V4-V6 denote the respective values of the genes, V7 is the gene ID (e.g., Ensembl ID), V8 is the gene name, and V9 is the distance of peak to nearest gene.
#' 
#' @param organism Object name assigned from readGFF() command.
#' @param start Lower bound of upstream extension.
#' @param end Upper bound of upstream extension.
#'
#' @return A data.table of unique genes located under peaks between two upstream extension levels.
#'
#' @examples
#' library(rtracklayer)
#' rat <- readGFF("ftp://ftp.ensembl.org/pub/release-84/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.84.gtf.gz")
#' fpath <- system.file("extdata", "somepeaksfile.txt", package="geneXtendeR")
#' peaksInput(fpath)
#' distinct(rat, 2000, 3000)
#'
#' @useDynLib geneXtendeR, .registration = TRUE
#' 
#' @import data.table
#'
#' @export
distinct <- function(organism, start, end) {
 if(!file.exists("peaks.txt")){
   stop("Please run peaksInput() function first!  See ?peaksInput for more information")
 } else {
        oopts = options(warn=-1)
        on.exit(options(oopts))
        run2 <- function(f1, f2, peakslist) {
            .C("extractpeaks", f1, f2, peakslist)[[3]]
        }
        
    		sapply(c(start, end), .geneXtender, organism, FALSE)
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
        col.names <- c("Chromosome", "Peak-Start", "Peak-End", "Gene-Chr", "Gene-Start", "Gene-End", "Gene-ID", "Gene-Name", "Distance")
        DT <- data.table::as.data.table(do.call("rbind", strsplit(finalList, split = "\t")))
        names(DT) <- col.names
        return(DT)
	}
}