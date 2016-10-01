#' Produces bar charts.
#'
#' Makes bar graphs showing the number of genes under peaks at various upstream extension levels.
#'
#' @param organism Object name assigned from readGFF() command.
#' @param start Lower bound of upstream extension.
#' @param end Upper bound of upstream extension.
#' @param by Interval between consecutive extensions.
#'
#' @return Creates bar charts.
#'
#' @examples
#' rat <- readGFF("ftp://ftp.ensembl.org/pub/release-84/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.84.gtf.gz")
#' fpath <- system.file("extdata", "somepeaksfile.txt", package="geneXtendeR")
#' peaksInput(fpath)
#' barChart(rat, 1000, 3000, 100)
#'
#' @importFrom dplyr filter select arrange
#' @import utils
#' @import graphics
#'
#' @useDynLib geneXtendeR, .registration = TRUE
#'
#'
#' @export
barChart <- function(organism, start, end, by) {
 if(!file.exists("peaks.txt")){
   message("Please run peaksInput() function first!  See ?peaksInput for more information")
 } else {
	oopts = options(warn = -1)
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
        levels(geneXtender.file$seqid) <- gsub("X", "100", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("chrX", "100", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("Y", "200", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("chrY", "200", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("MT", "300", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("chrMT", "300", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("MtDNA", "300", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("chrMtDNA", "300", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("M", "300", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("chrM", "300", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("Mito", "300", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("chrMito", "300", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("I", "1", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("chr1", "1", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("II", "2", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("chr2", "2", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("III", "3", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("chr3", "3", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("IV", "4", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("chr4", "4", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("V", "5", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("chr5", "5", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("VI", "6", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("chr6", "6", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("VII", "7", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("chr7", "7", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("VIII", "8", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("chr8", "8", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("IX", "9", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("chr9", "9", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("XI", "11", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("chr11", "11", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("XII", "12", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("chr12", "12", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("XIII", "13", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("chr13", "13", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("XIV", "14", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("chr14", "14", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("XV", "15", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("chr15", "15", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("XVI", "16", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("chr16", "16", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("chr17", "17", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("chr18", "18", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("chr19", "19", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("chr20", "20", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("chr21", "21", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("chr22", "22", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("chr23", "23", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("chr24", "24", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("chr25", "25", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("chr26", "26", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("chr27", "27", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("chr28", "28", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("chr29", "29", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("chr30", "30", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("chr31", "31", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("chr32", "32", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("chr33", "33", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("chr34", "34", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("chr35", "35", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("chr36", "36", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("chr37", "37", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("chr38", "38", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("chr39", "39", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("chr40", "40", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("chr41", "41", levels(geneXtender.file$seqid))
        levels(geneXtender.file$seqid) <- gsub("chr42", "42", levels(geneXtender.file$seqid))
        geneXtender.file$seqid = as.numeric(as.character(geneXtender.file$seqid))
        geneXtender.sorted <- dplyr::arrange(geneXtender.file, as.numeric(seqid), start)
        write.table(geneXtender.sorted, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE, sprintf("geneXtender_gtf_%s.bed", upstream))
    }
    
    
	sapply(seq(start, end, by), geneXtender)
	xlabs <- as.character(seq(start, end, by))
	xDeltas <- vapply(seq_along(xlabs)[-1], function(i) paste(xlabs[(i-1):i], collapse = "-"), "")
	gxFiles <- sprintf("geneXtender_gtf_%s.bed", seq(start, end, by))

    
    run <- function(f1, f2, num) {
        .C("extractnumber", f1, f2, num)[[3]]
    }
    
    num = 0
    numvec <- numeric()
	for (x in gxFiles) {
        numbers <- run(f1 = "peaks.txt", f2 = x, as.integer(num))
        numvec <- append(numvec, numbers)
	}
	
	barplot(numvec, names.arg = xlabs)
	
	on.exit(options(oopts))
	
	}
}
