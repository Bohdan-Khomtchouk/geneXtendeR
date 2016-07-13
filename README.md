[![Travis-CI Build Status](https://travis-ci.org/Bohdan-Khomtchouk/geneXtendeR.svg?branch=master)](https://travis-ci.org/Bohdan-Khomtchouk/geneXtendeR)
[![AUR](https://img.shields.io/aur/license/yaourt.svg?maxAge=2592000)]()
[![GitHub version](https://badge.fury.io/gh/Bohdan-Khomtchouk%2FgeneXtendeR.svg)](https://badge.fury.io/gh/Bohdan-Khomtchouk%2FgeneXtendeR)
# geneXtendeR

`geneXtendeR` is an R package for histone modification ChIP-seq analysis.  It is designed to optimally annotate a histone modification ChIP-seq peak input file with functionally important genomic features (e.g., genes associated with peaks) based on optimization calculations.  `geneXtendeR` optimally extends the boundaries of every gene in a genome by some genomic distance (in DNA base pairs) for the purpose of flexibly incorporating cis-regulatory elements, such as enhancers and promoters, as well as downstream elements that are important to the function of the gene (relative to an epigenetic histone modification ChIP-seq dataset). `geneXtendeR` computes optimal gene extensions tailored to the broadness of the specific epigenetic mark (e.g., H3K9me1, H3K27me3), as determined by a user-supplied ChIP-seq peak input file. As such, `geneXtendeR` maximizes the signal-to-noise ratio of locating genes closest to and directly under peaks. By performing a computational expansion of this nature, ChIP-seq reads that would initially not map strictly to a specific gene can now be optimally mapped to the regulatory regions of the gene, thereby implicating the gene as a potential candidate, and thereby making the ChIP-seq experiment more successful. Such an approach becomes particularly important when working with epigenetic histone modifications that have inherently broad peaks.

### Workflow

Sample biological workflow using `geneXtendeR` in combination with existing statistical software to analyze peak significance:

<img width="772" alt="workflow" src="https://cloud.githubusercontent.com/assets/9893806/15252937/0a043428-18fe-11e6-8951-2dbdeaa97fb8.png">

Subsequent gene ontology or network analysis may be conducted on genes associated with statistically significant peaks.

### Installation instructions

You can install the current GitHub version using the [devtools](https://github.com/hadley/devtools) package and the following two commands in R (since `gtf` is a big data package, it may take a few minutes to install, depending on the speed of your Internet connection):
```R
if (!require("devtools")) install.packages("devtools")
devtools::install_github("Bohdan-Khomtchouk/gtf")
devtools::install_github("Bohdan-Khomtchouk/geneXtendeR")
```
And then load the packages:
```R
library(gtf)
library(geneXtendeR)
```

Also, as specified in the `SystemRequirements` field of the `DESCRIPTION` file, `geneXtendeR` requires the installation of an external program called `bedtools`.  This program must be pre-installed on your computer prior to using `geneXtendeR`.  Detailed installation instructions can be found here: http://bedtools.readthedocs.io/en/latest/content/installation.html.  You may follow this recipe:  

```shell
$ wget https://github.com/arq5x/bedtools2/releases/download/v2.25.0/bedtools-2.25.0.tar.gz
$ tar -zxvf bedtools-2.25.0.tar.gz
$ cd bedtools2
$ make
$ cd bin
$ cp * /usr/local/bin
```

After `bedtools` has been installed, the `geneXtendeR` package is fully configured and setup for use.  You may check that `bedtools` has been successfully installed via:

```shell
$ bedtools
bedtools: flexible tools for genome arithmetic and DNA sequence analysis.
usage:    bedtools <subcommand> [options]
```

### Usage instructions
See `geneXtendeR.pdf` in `/vignettes` directory.

