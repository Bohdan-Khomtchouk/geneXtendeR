<div align="center">

# geneXtendeR

<img src="https://cloud.githubusercontent.com/assets/9893806/16898879/d14647c6-4bba-11e6-93d5-90bc802ac8e9.png">

##### R/Bioconductor package for optimized functional annotation of ChIP-seq data

![AUR](https://img.shields.io/aur/license/yaourt.svg)
[![Bioconductor Time](http://bioconductor.org/shields/years-in-bioc/geneXtendeR.svg)](http://bioconductor.org/packages/release/bioc/html/geneXtendeR.html "Bioconductor status")
<!-- [![Bioconductor Downloads](http://bioconductor.org/shields/downloads/geneXtendeR.svg)](https://bioconductor.org/packages/stats/bioc/geneXtendeR/ "Percentile downloads") -->

</div>

To see `geneXtendeR` in action (and why you may want to use it as part of your next ChIP-seq workflow), see Section 1 of the [vignette](https://github.com/Bohdan-Khomtchouk/geneXtendeR/tree/master/inst/doc). 

## Github installation instructions (for latest updated version)

You can install the current GitHub version using the [devtools](https://github.com/hadley/devtools) package:

```R
if (!require("devtools")) install.packages("devtools")
devtools::install_github("Bohdan-Khomtchouk/geneXtendeR")
```
And then load the package:

```R
library(geneXtendeR)
```

## Bioconductor landing page (updated biannually -- so likely out-of-date)

`geneXtendeR` can be found at: https://bioconductor.org/packages/geneXtendeR/.  This URL will redirect to the release landing page of the package and is the URL that should be used (in publications, etc.) to refer to `geneXtendeR`.  You can also refer specifically to the development, release, or specific numbered version of Bioconductor:

https://bioconductor.org/packages/devel/geneXtendeR/

https://bioconductor.org/packages/release/geneXtendeR/

https://bioconductor.org/packages/3.6/geneXtendeR/

## Bioconductor installation instructions

```R
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("geneXtendeR")
```

## Citation

If you are using `geneXtendeR` in your work, please cite the paper (http://dx.doi.org/10.1101/082347) accordingly.  **Note:** this paper is not up-to-date with latest developments.  Please stay tuned for updated preprint in 2018.  
