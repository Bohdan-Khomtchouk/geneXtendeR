<div align="center">

# geneXtendeR

<img src="https://cloud.githubusercontent.com/assets/9893806/16898879/d14647c6-4bba-11e6-93d5-90bc802ac8e9.png">

##### R/Bioconductor package for optimized functional annotation of ChIP-seq data

[![Travis-CI Build Status](https://travis-ci.org/Bohdan-Khomtchouk/geneXtendeR.svg?branch=master)](https://travis-ci.org/Bohdan-Khomtchouk/geneXtendeR)
[![AUR](https://img.shields.io/aur/license/yaourt.svg?maxAge=2592000)]()
[![Bioconductor Time](http://bioconductor.org/shields/years-in-bioc/geneXtendeR.svg)](http://bioconductor.org/packages/release/bioc/html/geneXtendeR.html "Bioconductor status")
[![Bioconductor Downloads](http://bioconductor.org/shields/downloads/geneXtendeR.svg)](http://bioconductor.org/packages/stats/bioc/geneXtendeR.html "Percentile downloads")
[![Bioconductor Commits](http://bioconductor.org/shields/commits/bioc/geneXtendeR.svg)](http://bioconductor.org/packages/release/bioc/html/geneXtendeR.html#svn_source "svn commits")
[![Support posts](http://bioconductor.org/shields/posts/geneXtendeR.svg)](https://support.bioconductor.org/t/geneXtendeR/ "Bioconductor support posts")

</div>

`geneXtendeR` optimizes the functional annotation of ChIP-seq peaks using fast iterative peak-coordinate/GTF alignment algorithms.  Since different ChIP-seq peak callers produce different differentially enriched peaks with a large variance in peak length distribution and total peak count, annotating peak lists with their nearest genes can be a noisy process.  As such, the goal of `geneXtendeR` is to robustly link differentially enriched peaks with their respective genes, thereby aiding experimental follow-up and validation in designing primers for a set of prospective gene candidates during qPCR.

### Landing page

`geneXtendeR` can be found at: https://bioconductor.org/packages/geneXtendeR/.  This URL will redirect to the release landing page of the package and is the URL that should be used (in publications, etc.) to refer to `geneXtendeR`.  You can also refer specifically to the development, release, or specific numbered version of Bioconductor:

https://bioconductor.org/packages/devel/geneXtendeR/

https://bioconductor.org/packages/release/geneXtendeR/

https://bioconductor.org/packages/3.6/geneXtendeR/

### Vignette

Please read the `geneXtendeR` documentation (https://bioconductor.org/packages/devel/bioc/vignettes/geneXtendeR/inst/doc/geneXtendeR.pdf), which includes software workflows and comprehensive analyses to get you started.

### Github installation instructions

You can install the current GitHub version using the [devtools](https://github.com/hadley/devtools) package:

```R
if (!require("devtools")) install.packages("devtools")
devtools::install_github("Bohdan-Khomtchouk/geneXtendeR")
```
And then load the package:

```R
library(geneXtendeR)
```

### Bioconductor installation instructions

```R
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("geneXtendeR")
```

After installation, you should always update to the latest available version of `geneXtendeR` by typing `biocLite()`:

```R
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite()                  ## R version 3.0 or later
```

### Citation

If you are using `geneXtendeR` in your work, please cite the paper (http://dx.doi.org/10.1101/082347) accordingly.
