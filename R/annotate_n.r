#' Annotate peaks file.
#'
#' Annotate a user's peaks file (which has been preprocessed with the peaksInput() command) with gene information based on optimally chosen geneXtendeR upstream extension file.  This command requires a preprocessed "peaks.txt" file (generated using peaksInput()) to be present in the user's working directory, otherwise the user is prompted to rerun the peaksInput() command in order to regenerate it.
#' 
#' @param organism Object name assigned from readGFF() command.
#' @param extension Desired upstream extension.
#' @param n Number of Gene's closest away from the peak
#'
#' @return The gene coordinates are extended by `extension` at the 5-prime end, and by 500 bp at the 3-prime end.  The peaks file is then overlayed on these new gene coordinates, producing a file of peaks annotated with gene ID, gene name, and gene-to-peak genomic distance (in bp).  Distance is calculated between 5-prime end of gene and 3-prime end of peak. File named "annotated_'extension'_'n'.txt".
#' @return A data.table formatted version of the annotated file for checking or further calculations.
#'
#' @examples
#' rat <- readGFF("ftp://ftp.ensembl.org/pub/release-84/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.84.gtf.gz")
#' fpath <- system.file("extdata", "somepeaksfile.txt", package="geneXtendeR")
#' peaksInput(fpath)
#' annotate_n(rat, 2500, n=3)
#'
#' @useDynLib geneXtendeR, .registration = TRUE
#'
#' @export
annotate_n <- function(organism, extension, n=2) {
  if(!file.exists("peaks.txt")){
    message("Please run peaksInput() function first!  See ?peaksInput for more information")
  } else if(n < 1) {
    message("Please run a valid distance of gene's away with n > 0")
  } else {
    oopts = options(warn=-1)
    on.exit(options(oopts))
    
    genes <- .geneXtender(extension, organism, TRUE)
    
    run3 <- function(f1, f2, peakslist) {
      .C("annotate", f1, f2, peakslist)[[3]]
    }
    
    linelen <- "                                                                                                    "
    peaksArray<-rep(linelen,500000)
    onegxFile <- sprintf("geneXtender_gtf_%s.bed", extension)
    onecmd2 <- run3(f1 = "peaks.txt", f2 = onegxFile, peakslist = peaksArray) 
    onecmd3 <- onecmd2[onecmd2 != linelen]
    
    rdt <- data.table::transpose(data.table::as.data.table(strsplit(onecmd3, "\t"))) #Formatting
    colnames(rdt) <- c("Chromosome", "Peak-Start", "Peak-End", "rm", "Gene-Start", "Gene-End", "Gene-ID", "Gene-Name", "Distance-of-Gene-to-Nearest-Peak")
    rdt <- rdt[,rm:= NULL]
    for (j in c("Chromosome", "Peak-Start", "Peak-End", "Gene-Start", "Gene-End", "Distance-of-Gene-to-Nearest-Peak")) data.table::set(rdt, j=j, value = as.numeric(rdt[[j]]))

    
    #indexing
    indices <- sapply(rdt$`Gene-ID`, function(x) which(x == genes$gene_id))
    
    #error_checking (MORETODO)
    indices <- ifelse(indices < (n), n, indices)
    indices <- ifelse(indices > (length(genes$seqid) - (n-1)), length(genes$seqid) - (n-1), indices)
    rdt[,'bot' := (indices - (n-1))]
    rdt[,'top' := (indices + (n-1))]
    
    #THROW extras
    rdt[, c("Gene-Start", "Gene-End", "Gene-ID", "Gene-Name", "Distance-of-Gene-to-Nearest-Peak") := NULL]
    
    #Grabs and adds
    rdt[, ..I := .I]
    some_fun <- function(dtx) {
      g <- genes[dtx$bot:dtx$top]
      g[,names(dtx) := dtx]
      return(g)
    }
    rdt <- rdt[, some_fun(.SD), by = ..I] #applies function by row
    rdt[, c('top', 'bot') := NULL] #Remove top and bot info
    
    es <- rdt$`Peak-End` - rdt$start
    se <- rdt$`Peak-Start` - rdt$end
    
    rdt[, "Minimum-Distance-to-Gene" := ifelse(between(0, pmin(es, se), pmax(es, se)), 0 ,pmin(abs(es), abs(se)))]
    
    setorderv(rdt, c("..I", "Minimum-Distance-to-Gene"))
    rdt[,rank := rep(1:(n*2-1), (.N/(n*2-1)))]
    rdt[rank <= n]
    
    
    data.table::setnames(rdt, c("start", "end", "gene_id", "gene_name", "..I"), c("Gene-Start", "Gene-End", "Gene-ID", "Gene-Name", "Peak-Num"))
    data.table::setcolorder(rdt, c("Peak-Num", "Chromosome", "Peak-Start", "Peak-End", "Gene-Start", "Gene-End", "Gene-ID", "Gene-Name", "rank", "Minimum-Distance-to-Gene", "seqid"))
    
    #TODO more error handling (for now)
    err <- which(rdt$Chromosome != rdt$seqid)
    # if (length(err) > 0) {
    #   message(sprintf("Warning! There were %s peaks that could not find the next gene on their respective chromosomes %s genes away", length(err), n))
    #   #for(col in c("seqid", "start", "end", "gene_id", "gene_name")) set(export_table, i=err, j=col, value=NA) #cleanup
    # }
    rdt <- rdt[rank<= n]
    data.table::fwrite(
      rdt,
      file = sprintf("annotated_%s_%s.txt", extension, n),
      sep = "\t",
      row.names = FALSE,
      col.names = TRUE,
      quote = FALSE
    )
    return(rdt)
  }
}