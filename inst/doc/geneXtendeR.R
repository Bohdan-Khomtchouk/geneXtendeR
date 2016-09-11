### R code from vignette source 'geneXtendeR.Rnw'

###################################################
### code chunk number 1: <style-Sweave
###################################################
BiocStyle::latex()


###################################################
### code chunk number 2: geneXtendeR.Rnw:33-37
###################################################
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("geneXtendeR")
library(geneXtendeR)


###################################################
### code chunk number 3: geneXtendeR.Rnw:42-44 (eval = FALSE)
###################################################
## rat <- readGFF("ftp://ftp.ensembl.org/pub/release-84/gtf/
##                       rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.84.chr.gtf.gz")


###################################################
### code chunk number 4: geneXtendeR.Rnw:51-52
###################################################
peaksInput("somepeaksfile.txt")


###################################################
### code chunk number 5: geneXtendeR.Rnw:61-62
###################################################
barChart(rat, 0, 10000, 500)


###################################################
### code chunk number 6: geneXtendeR.Rnw:73-74
###################################################
linePlot(rat, 0, 10000, 500)


###################################################
### code chunk number 7: geneXtendeR.Rnw:83-84
###################################################
linePlot(rat, 2000, 3000, 100)


###################################################
### code chunk number 8: geneXtendeR.Rnw:94-95 (eval = FALSE)
###################################################
## cat(distinct(rat, 2300, 2400))


###################################################
### code chunk number 9: geneXtendeR.Rnw:112-113
###################################################
linePlot(rat, 7000, 8500, 100)


###################################################
### code chunk number 10: geneXtendeR.Rnw:123-124 (eval = FALSE)
###################################################
## cat(distinct(rat, 7800, 7900))


###################################################
### code chunk number 11: geneXtendeR.Rnw:133-134
###################################################
linePlot(rat, 0, 10000, 500)


###################################################
### code chunk number 12: geneXtendeR.Rnw:141-142
###################################################
annotate(rat, 3300)


