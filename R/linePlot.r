utils::globalVariables(c("gxFiles", "xDeltas"))
#' Produces line plots.
#'
#' Makes differential line plots showing the differences in the number of genes under peaks at consecutive upstream extension levels.
#'
#' @return Creates differential line plots.
#'
#' @examples
#' \dontrun{generate(rat, 1000, 3000, 100)}
#' \dontrun{linePlot()}  #Explained in vignette
#'
#' @export
linePlot <- function() {
		options(warn = -1)
		for (x in gxFiles) {
			cmds <- sprintf("bedtools closest -d -a peaks.txt -b %s", gxFiles)
		}
		countZeros <- function(y) {
			bedtools.closest.output <- system(cmds[y], intern = TRUE, ignore.stderr = TRUE)
			length(grep("\\<0\\>", unlist(strsplit(bedtools.closest.output, "\t"))))
		}
		z <- sapply(1:length(cmds), countZeros)
		diffs <- sapply(1:length(z), function(i) {z[i+1] - z[i]})
		differences <- diffs[!is.na(diffs)]
		par(mar = c(7.1,4.1,2.1,2.1))
		plot(differences, type = "o", col = "blue", xaxt = "n", ylab = "differences", xlab = "")
		axis(1, at = 1:length(differences), labels = xDeltas, las = 3)
		mtext(side = 1, "Genomic region (bp)", line = 6)
	}
