#' Looks up specific gene and closest peaks
#'
#' Looks up closest peak to a specified gene on the peaks file (which has been preprocessed with the peaksInput() command) based on the latest .bed file accessed or for a specified extension.  This command requires a preprocessed "peaks.txt" file (generated using peaksInput()) to be present in the user's working directory, otherwise the user is prompted to rerun the peaksInput() command in order to regenerate it.
#' 
#' @param organism Object name assigned from readGFF() command.
#' @param gene_name Gene names or gene ids specified by user in a string form. 
#' @param n Number of closest peaks located to `gene_name` on any given chromosome to be found. Default = 2
#' @param extension Desired upstream extension. If left unspecified, the latest geneXtender bed file will be chosen. If no extension is specified and no bed file can be found, a default extension of 500 is selected.
#' @param cutoff Optional arg to specify max bp distance to search around `gene_name`. Default = Inf
#'
#' @return A data.table with all peaks located within `n` peaks and `cutoff` bp distance away on every chromosome that contains `gene_name`.
#'
#' @examples
#' library(rtracklayer)
#' rat <- readGFF("ftp://ftp.ensembl.org/pub/release-84/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.84.gtf.gz")
#' fpath <- system.file("extdata", "somepeaksfile.txt", package="geneXtendeR")
#' closest <- gene_lookup(rat, c("Vom2r3", "Vom2r5"), n = 7, extension = 1000)
#' closest
#'
#' @useDynLib geneXtendeR, .registration = TRUE
#'
#' @export
gene_lookup <- function(organism, gene_name, n=2, extension = NA, cutoff = Inf) {
  if(!file.exists("peaks.txt")){
    message("Please run peaksInput() function first!  See ?peaksInput for more information")
  } else if (cutoff < 0 | n < 1) {
    message("Please select a cutoff >= 0 or a value of n >= 1")
  } else {
    genes <- NULL
    
    if (!is.na(extension)) {
      file <- sprintf("geneXtender_gtf_%i.bed", extension)
      # if (file.exists(file)) { #Only useful assuming that the organisms are the same. unable to check so leave commented
      #   genes <- data.table::fread(file)
      # }
      # else {
      genes <- .geneXtender(extension, organism, TRUE)
      # }
    } else {
      files <- grep("geneXtender_gtf_*", list.files(), value=TRUE)
      if (length(files) == 0) { #If no bed file exists, generate one with default upstream of 500
        genes <- .geneXtender(500, organism, TRUE)
        files <- grep("geneXtender_gtf_*", list.files(), value=TRUE)
        files_check <- c("geneXtender_gtf_500.bed")
      }
      
      file_info <- file.info(files)
      file_info <- file_info[order(-as.numeric(file_info$atime)),]
      recent_extension <- row.names(file_info)[1]
      genes <- data.table::fread(recent_extension)
    }
    colnames(genes) <- c("seqid", "start", "end", "gene_id", "gene_name_id")
    
    genez <- genes[tolower(gene_name_id) %in% tolower(gene_name) | tolower(gene_id) %in% tolower(gene_name),]
    
    peaks <- data.table::fread("peaks.txt")
    colnames(peaks) <- c("Chromosome", "Peak-Start", "Peak-End")
    
    internal_find <- function(dtc) {
      peak <- subset(peaks, Chromosome == dtc$seqid)
      
      peak$es <- dtc$start - peak$`Peak-End`
      peak$se <- dtc$end - peak$`Peak-Start`
      
      peak$distance <- ifelse(data.table::between(0, pmin(peak$es, peak$se), pmax(peak$es,peak$se)), 0, pmin(abs(peak$es), abs(peak$se)))
      
      peak <- data.table::setorderv(peak, "distance")
      
      zeroes <- peak[, distance == 0]
      
      vals <- peak[order(distance)]
      vals <- na.omit(vals[1:n,])
      
      
      #ERROR CHECKING
      if(nrow(peak) < n) {
        warning(sprintf("Warning! There were more peaks requested on Chromosome on %i than exist on the peaks.txt file.
                        Only %i peaks could be found on Chromosome %i and were returned.", dtc$seqid, nrow(peak), dtc$seqid))
      }
      return(cbind(vals, dtc))
    }
    
    genez[, ..I := .I]
    report <- genez[, internal_find(.SD), by = ..I]
    report <- report[distance <= cutoff]
    
    report[, c('..I', 'es', 'se', 'seqid', 'gene_id') := NULL]
    data.table::setnames(report, c("Chromosome", "Peak-Start", "Peak-End", "Distance-to-Gene", "Gene-Start", "Gene-End", "Gene"))
    
    
    return(report[])
  }
}
