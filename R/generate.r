utils::globalVariables(c("V3", "V9", "V1", "V4", "V5", "V7", "gene_id", "gene_name", "xlabs", "xDeltas", "gxFiles"))
#' Generate geneXtendeR files.
#'
#' Creates organism-specific customized gene transfer format files at user-designated extensions.  Only upstream extensions can be specified.  Downstream extensions are always set as default to 500 bp (to capture any downstream gene elements).
#'
#' @param organism Name of organism (refer to gtf R package).
#' @param start Lower bound of upstream extension.
#' @param end Upper bound of upstream extension.
#' @param by Interval between consecutive extensions.
#'
#' @return Generates exactly ([(start-end)/by] + 1) modified gtf files and loads in global variables for subsequent analysis.
#'
#' @examples
#' \dontrun{generate(rat, 2000, 10000, 500)}
#'
#' @import data.table
#' @import gtf
#' 
#' @export
generate <- function(organism, start, end, by) {

  options(warn = -1)
  geneXtender <- function(upstream) {
    messy2 <- dplyr::filter(organism, V3 == "gene")
    messy3 <- tidyr::extract(messy2, V9, c('gene_id', 'gene_name'), 'gene_id .(\\S+).;.*gene_name .(\\S+).;')
    neat <- dplyr::select(messy3, V1, V4, V5, V7, gene_id, gene_name)
    pos_exons <- dplyr::filter(neat, V7 == "+")
    neg_exons <- dplyr::filter(neat, V7 == "-")
    pos_exons$V4 = pos_exons$V4 - upstream
    pos_exons$V4[pos_exons$V4 < 0] <- 1
    pos_exons$V5 = pos_exons$V5 + 500
    neg_exons$V4 = neg_exons$V4 - 500
    neg_exons$V4[neg_exons$V4 < 0] <- 1
    neg_exons$V5 = neg_exons$V5 + upstream
    merged_exons <- rbind(pos_exons, neg_exons)
    geneXtender.file <- dplyr::select(merged_exons, V1, V4, V5, gene_id, gene_name)
    setDT(geneXtender.file)
    geneXtender.file[V1 == "X", V1 := 100]
    geneXtender.file[V1 == "Y", V1 := 200]
    geneXtender.file[V1 == "MT", V1 := 300]
    geneXtender.file[V1 == "MtDNA", V1 := 300]
    geneXtender.file[V1 == "M", V1 := 300]
    geneXtender.file[V1 == "Mito", V1 := 300]
    geneXtender.file[V1 == "I", V1 := 1]
    geneXtender.file[V1 == "II", V1 := 2]
    geneXtender.file[V1 == "III", V1 := 3]
    geneXtender.file[V1 == "IV", V1 := 4]
    geneXtender.file[V1 == "V", V1 := 5]
    geneXtender.file[V1 == "VI", V1 := 6]
    geneXtender.file[V1 == "VII", V1 := 7]
    geneXtender.file[V1 == "VIII", V1 := 8]
    geneXtender.file[V1 == "IX", V1 := 9]
    geneXtender.file[V1 == "XI", V1 := 11]
    geneXtender.file[V1 == "XII", V1 := 12]
    geneXtender.file[V1 == "XIII", V1 := 13]
    geneXtender.file[V1 == "XIV", V1 := 14]
    geneXtender.file[V1 == "XV", V1 := 15]
    geneXtender.file[V1 == "XVI", V1 := 16]
    geneXtender.sorted <- dplyr::arrange(geneXtender.file, as.numeric(V1), V4)
    write.table(geneXtender.sorted, quote = FALSE, sep = ",", row.names = FALSE, col.names = FALSE, sprintf("geneXtender_gtf_%s.txt", upstream))
  }

	sapply(seq(start, end, by), geneXtender)
	assign("xlabs", as.character(seq(start, end, by)), envir = .GlobalEnv)
	assign("xDeltas", vapply(seq_along(xlabs)[-1], function(i) paste(xlabs[(i-1):i], collapse = "-"), ""), envir = .GlobalEnv)
	assign("gxFiles", sprintf("geneXtender_gtf_%s.txt", seq(start, end, by)), envir = .GlobalEnv)

}
