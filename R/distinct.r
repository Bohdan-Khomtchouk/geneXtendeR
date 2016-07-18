utils::globalVariables(c("type", "seqid", "start", "end", "strand", "gene_id", "gene_name", "cmd1", "cmd2", "cmd2.bedtools.closest.output.zeros.DT"))
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
#' @return A table of unique genes located under peaks between two upstream extension levels.
#'
#' @examples
#' \dontrun{rat <- readGFF("ftp://ftp.ensembl.org/pub/release-84/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.84.gtf.gz")}
#' \dontrun{generate(rat, 1000, 3000, 500)}
#' \dontrun{distinct(rat, 2000, 3000)}
#'
#' @export
distinct <- function(organism, start, end) {
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

		sapply(c(start, end), geneXtender)
		twogxFiles <- sprintf("geneXtender_gtf_%s.bed", c(start, end))

		for (m in twogxFiles) {
			twocmds <- sprintf("bedtools closest -d -a peaks.txt -b %s", twogxFiles)
			assign("cmd1", twocmds[1], envir = .GlobalEnv)
			assign("cmd2", twocmds[2], envir = .GlobalEnv)
		}

		cmd1.bedtools.closest.output <- system(cmd1, intern = TRUE, ignore.stderr = TRUE)
		cmd2.bedtools.closest.output <- system(cmd2, intern = TRUE, ignore.stderr = TRUE)
		cmd1.bedtools.closest.output.zeros <- subset(cmd1.bedtools.closest.output, grepl("^.+(\t0)$", cmd1.bedtools.closest.output))
		cmd2.bedtools.closest.output.zeros <- subset(cmd2.bedtools.closest.output, grepl("^.+(\t0)$", cmd2.bedtools.closest.output))
		cmd1.bedtools.closest.output.zeros.DT <- as.data.table(do.call("rbind", strsplit(cmd1.bedtools.closest.output.zeros, split = "\t")))
		assign("cmd2.bedtools.closest.output.zeros.DT", as.data.table(do.call("rbind", strsplit(cmd2.bedtools.closest.output.zeros, split = "\t"))), envir = .GlobalEnv)
		cmd2.bedtools.closest.output.zeros.DT[!cmd1.bedtools.closest.output.zeros.DT, on = paste0("V", 1:3)]

}