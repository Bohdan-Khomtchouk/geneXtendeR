utils::globalVariables(c("type", "seqid", "start", "end", "strand", "gene_id", "gene_name", "xlabs", "xDeltas", "gxFiles"))
#' Generate geneXtendeR files.
#'
#' Creates organism-specific customized gene transfer format files at user-designated extensions.  Only upstream extensions can be specified.  Downstream extensions are always set as default to 500 bp (to capture any downstream gene elements).
#'
#' @param organism Object name assigned from readGFF() command.
#' @param start Lower bound of upstream extension.
#' @param end Upper bound of upstream extension.
#' @param by Interval between consecutive extensions.
#'
#' @return Generates exactly ([(start-end)/by] + 1) modified gtf files and loads in global variables for subsequent analysis.
#'
#' @examples
#' \dontrun{human <- readGFF("ftp://ftp.ensembl.org/pub/release-84/gtf/homo_sapiens/Homo_sapiens.GRCh38.84.gtf.gz")}
#' \dontrun{generate(human, 2000, 10000, 500)}
#'
#' @import rtracklayer
#' @import data.table
#' 
#' @export
generate <- function(organism, start, end, by) {

  options(warn = -1)
  geneXtender <- function(upstream) {
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
    setDT(geneXtender.file)
    geneXtender.file[seqid == "X", seqid := 100]
    geneXtender.file[seqid == "Y", seqid := 200]
    geneXtender.file[seqid == "MT", seqid := 300]
    geneXtender.file[seqid == "MtDNA", seqid := 300]
    geneXtender.file[seqid == "M", seqid := 300]
    geneXtender.file[seqid == "Mito", seqid := 300]
    geneXtender.file[seqid == "I", seqid := 1]
    geneXtender.file[seqid == "II", seqid := 2]
    geneXtender.file[seqid == "III", seqid := 3]
    geneXtender.file[seqid == "IV", seqid := 4]
    geneXtender.file[seqid == "V", seqid := 5]
    geneXtender.file[seqid == "VI", seqid := 6]
    geneXtender.file[seqid == "VII", seqid := 7]
    geneXtender.file[seqid == "VIII", seqid := 8]
    geneXtender.file[seqid == "IX", seqid := 9]
    geneXtender.file[seqid == "XI", seqid := 11]
    geneXtender.file[seqid == "XII", seqid := 12]
    geneXtender.file[seqid == "XIII", seqid := 13]
    geneXtender.file[seqid == "XIV", seqid := 14]
    geneXtender.file[seqid == "XV", seqid := 15]
    geneXtender.file[seqid == "XVI", seqid := 16]
    geneXtender.sorted <- dplyr::arrange(geneXtender.file, as.numeric(seqid), start)
    write.table(geneXtender.sorted, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE, sprintf("geneXtender_gtf_%s.bed", upstream))
  }

	sapply(seq(start, end, by), geneXtender)
	assign("xlabs", as.character(seq(start, end, by)), envir = .GlobalEnv)
	assign("xDeltas", vapply(seq_along(xlabs)[-1], function(i) paste(xlabs[(i-1):i], collapse = "-"), ""), envir = .GlobalEnv)
	assign("gxFiles", sprintf("geneXtender_gtf_%s.bed", seq(start, end, by)), envir = .GlobalEnv)

}
