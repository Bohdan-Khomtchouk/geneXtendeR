#' Function for formatting data to be called by C function.
#'
#' Formats genome-data based on upstream and organism for .C caller function. Not to be run by user.
#' 
#' @param upstream Extension upstream of each gene.
#' @param organism Object name assigned from readGFF() command.
#' @param ret Logical that designates if the output is written in R or only as a .bed file
#'
#' @return The gene coordinates extended by the designated upstream. In .bed file format.
#' @return A data.table formatted version of the data depending on the ret logical.
#'
#' @note Not to be run by the user
#' 
#' @keywords internal
#' 
#'
#' @useDynLib geneXtendeR, .registration = TRUE
.geneXtender <- function(upstream, organism, ret = FALSE) {
  oopts = options(warn=-1)
  on.exit(options(oopts))
  messy2 <- dplyr::filter(organism, type == "gene")
  neat <- dplyr::select(messy2, seqid, start, end, strand, gene_id, gene_name)
  pos_exons <- dplyr::filter(neat, strand == "+")
  neg_exons <- dplyr::filter(neat, strand == "-")
  pos_exons$start = pos_exons$start - upstream
  pos_exons$start[pos_exons$start < 0] <- 1
  pos_exons$end = pos_exons$end + 500
  neg_exons$start = neg_exons$start - 500
  neg_exons$start[neg_exons$start < 0] <- 1
  neg_exons$end = neg_exons$end + upstream
  merged_exons <- rbind(pos_exons, neg_exons)
  geneXtender.file <- dplyr::select(merged_exons, seqid, start, end, gene_id, gene_name)
  levels(geneXtender.file$seqid) <- gsub("chr", "", levels(geneXtender.file$seqid))
  levels(geneXtender.file$seqid) <- gsub("CHR", "", levels(geneXtender.file$seqid))
  levels(geneXtender.file$seqid) <- gsub("X", "100", levels(geneXtender.file$seqid))
  levels(geneXtender.file$seqid) <- gsub("Y", "200", levels(geneXtender.file$seqid))
  levels(geneXtender.file$seqid) <- gsub("MT", "300", levels(geneXtender.file$seqid))
  levels(geneXtender.file$seqid) <- gsub("MtDNA", "300", levels(geneXtender.file$seqid))
  levels(geneXtender.file$seqid) <- gsub("M", "300", levels(geneXtender.file$seqid))
  levels(geneXtender.file$seqid) <- gsub("Mito", "300", levels(geneXtender.file$seqid))
  levels(geneXtender.file$seqid) <- gsub("I", "1", levels(geneXtender.file$seqid))
  levels(geneXtender.file$seqid) <- gsub("II", "2", levels(geneXtender.file$seqid))
  levels(geneXtender.file$seqid) <- gsub("III", "3", levels(geneXtender.file$seqid))
  levels(geneXtender.file$seqid) <- gsub("IV", "4", levels(geneXtender.file$seqid))
  levels(geneXtender.file$seqid) <- gsub("V", "5", levels(geneXtender.file$seqid))
  levels(geneXtender.file$seqid) <- gsub("VI", "6", levels(geneXtender.file$seqid))
  levels(geneXtender.file$seqid) <- gsub("VII", "7", levels(geneXtender.file$seqid))
  levels(geneXtender.file$seqid) <- gsub("VIII", "8", levels(geneXtender.file$seqid))
  levels(geneXtender.file$seqid) <- gsub("IX", "9", levels(geneXtender.file$seqid))
  levels(geneXtender.file$seqid) <- gsub("XI", "11", levels(geneXtender.file$seqid))
  levels(geneXtender.file$seqid) <- gsub("XII", "12", levels(geneXtender.file$seqid))
  levels(geneXtender.file$seqid) <- gsub("XIII", "13", levels(geneXtender.file$seqid))
  levels(geneXtender.file$seqid) <- gsub("XIV", "14", levels(geneXtender.file$seqid))
  levels(geneXtender.file$seqid) <- gsub("XV", "15", levels(geneXtender.file$seqid))
  levels(geneXtender.file$seqid) <- gsub("XVI", "16", levels(geneXtender.file$seqid))
  geneXtender.file$seqid = as.numeric(as.character(geneXtender.file$seqid))
  geneXtender.sorted <- dplyr::arrange(geneXtender.file, as.numeric(seqid), start)
  write.table(geneXtender.sorted, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE, sprintf("geneXtender_gtf_%s.bed", upstream))
  if (ret == TRUE) {
    return(data.table::as.data.table(geneXtender.sorted))
  }
}