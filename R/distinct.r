utils::globalVariables(c("V3", "V9", "V1", "V4", "V5", "V7", "gene_id", "gene_name", "cmd1", "cmd2", "cmd2.bedtools.closest.output.zeros.DT"))
#' Finds unique genes under peaks.
#'
#' Determines what genes directly under peaks are actually unique between two different upstream extension levels.
#' 
#' V1-V3 denote the chromosome/start/end positions of the peaks, V4-V6 denote the respective values of the genes, V7 is the gene ID (e.g., Ensembl ID), V8 is the gene name, and V9 is the distance of peak to nearest gene.
#' 
#' @param organism Name of organism (refer to gtf R package).
#' @param start Lower bound of upstream extension.
#' @param end Upper bound of upstream extension.
#'
#' @return A table of unique genes located under peaks between two upstream extension levels.
#'
#' @examples
#' \dontrun{generate(rat, 1000, 3000, 500)}
#' \dontrun{distinct(rat, 2000, 3000)}
#'
#' @export
distinct <- function(organism, start, end) {
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

		sapply(c(start, end), geneXtender)
		twogxFiles <- sprintf("geneXtender_gtf_%s.txt", c(start, end))

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