#' Annotate peaks file.
#'
#' Annotate peaks file with gene information based on optimally chosen geneXtendeR upstream extension file.
#' 
#' @param organism Object name assigned from readGFF() command.
#' @param extension Desired upstream extension.
#'
#' @return A file of peaks annotated with gene ID, gene name, and gene-to-peak genomic distance (in bp).
#'
#' @examples
#' \dontrun{rat <- readGFF("ftp://ftp.ensembl.org/pub/release-84/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.84.gtf.gz")}
#' \dontrun{peaksInput("somepeaksfile.txt")}
#' \dontrun{annotate(rat, 2500)}
#'
#' @import data.table
#'
#' @export
annotate <- function(organism, extension) {
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
    			data.table::setDT(geneXtender.file)
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

		sapply(extension, geneXtender)
		onegxFile <- sprintf("geneXtender_gtf_%s.bed", extension)
		onecmd <- sprintf("bedtools closest -d -a peaks.txt -b %s", onegxFile)
		onecmd.bedtools.closest.output <- system(onecmd, intern = TRUE, ignore.stderr = TRUE)
		write.table(read.table(text = onecmd.bedtools.closest.output), 
		            file = sprintf("peaks_annotated_%s.txt", extension),
		            quote = FALSE, 
		            sep = "\t", 
		            row.names = FALSE,
		            col.names = c("Chromosome", "Peak Start", "Peak End", "Chromosome", "Gene Start", "Gene End", "Gene ID", "Gene Name", "Distance-of-Gene-to-Nearest-Peak")
		            )
}