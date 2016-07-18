utils::globalVariables(c("gxFiles", "xlabs"))
#' Produces bar charts.
#'
#' Makes bar graphs showing the number of genes under peaks at various upstream extension levels.
#'
#' @return Creates bar charts.
#'
#' @examples
#' \dontrun{mouse <- readGFF("ftp://ftp.ensembl.org/pub/release-84/gtf/mus_musculus/Mus_musculus.GRCm38.84.gtf.gz")}
#' \dontrun{generate(mouse, 1000, 3000, 100)}
#' \dontrun{barChart()}
#'
#' @export
barChart <- function() {
	options(warn = -1)
	for (x in gxFiles) {
		cmds <- sprintf("bedtools closest -d -a peaks.txt -b %s", gxFiles)
	}
	countZeros <- function(y) {
		bedtools.closest.output <- system(cmds[y], intern = TRUE, ignore.stderr = TRUE)
		length(grep("\\<0\\>", unlist(strsplit(bedtools.closest.output, "\t"))))
		}
	barPoints <- sapply(1:length(cmds), countZeros)
	barplot(barPoints, names.arg = xlabs)
	}
