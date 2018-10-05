#' Finds differential gene ontologies
#'
#' Determines gene ontology terms for each category (biological process (BP), cellular compartment (CC), molecular function (MF)) of genes-under-peaks that are unique between two different upstream extension levels.
#'
#'
#' @param organism Object name assigned from readGFF() command.
#' @param start Lower bound of upstream extension.
#' @param end Upper bound of upstream extension.
#' @param GOcategory Either BP, CC, or MF.
#' @param GOspecies Either org.Ag.eg.db (mosquito), org.Bt.eg.db (bovine), org.Ce.eg.db (worm), org.Cf.eg.db (canine), org.Dm.eg.db (fly), org.Dr.eg.db (zebrafish), org.Gg.eg.db (chicken), org.Hs.eg.db (human), org.Mm.eg.db (mouse), org.Mmu.eg.db (rhesus), org.Pt.eg.db (chimpanzee), org.Rn.eg.db (rat), org.Sc.sgd.db (yeast), org.Ss.eg.db (pig), or org.Xl.eg.db (frog).  
#'
#' @return A data frame of gene symbol, gene ontology ID, and gene ontology term for either a BP, CC, or MF category.  This data frame displays the annotations of all unique genes (i.e., genes that are located under peaks between two upstream extension levels) with their respective gene ontology information.
#'
#' @importFrom GO.db GO.db
#' 
#' @examples
#' library(rtracklayer)
#' library(org.Rn.eg.db)
#' rat <- readGFF("ftp://ftp.ensembl.org/pub/release-84/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.84.gtf.gz")
#' fpath <- system.file("extdata", "somepeaksfile.txt", package="geneXtendeR")
#' peaksInput(fpath)
#' diffGO(rat, 0, 500, BP, org.Rn.eg.db)
#'
#' @useDynLib geneXtendeR, .registration = TRUE
#'
#' @export

diffGO <- function(organism, start, end, GOcategory, GOspecies) {
  if(!file.exists("peaks.txt")){
    message("Please run peaksInput() function first!  See ?peaksInput for more information")
  } else {

    run2 <- function(f1, f2, peakslist) {
      .C("extractpeaks", f1, f2, peakslist)[[3]]
    }

    sapply(c(start, end), .geneXtender, rat, FALSE)
    twogxFiles <- sprintf("geneXtender_gtf_%s.bed", c(start, end))
    linelen <- ""  
    n <- 500000
    peaksArray <- rep(linelen, n)
    peaksArray2 <- rep(linelen, n)
    cmdtmp1 <- run2(f1 = "peaks.txt", f2 = twogxFiles[[1]], as.character(peaksArray))
    cmdtmp2 <- run2(f1 = "peaks.txt", f2 = twogxFiles[[2]], as.character(peaksArray2))
    cmd1 <- cmdtmp1[cmdtmp1 != linelen]
    cmd2 <- cmdtmp2[cmdtmp2 != linelen]
    m = regexec("^(?:[^\t]+\t){3}", cmd1)
    first3.cmd1 = unlist(regmatches(cmd1, m))
    m = regexec("^(?:[^\t]+\t){3}", cmd2)
    first3.cmd2 = unlist(regmatches(cmd2, m))
    finalList = cmd2[!(first3.cmd2 %in% first3.cmd1)]
    DT <- data.table::as.data.table(do.call("rbind", strsplit(finalList, split = "\t")))
    gene_names <- DT[[8]]

    if (deparse(substitute(GOcategory)) == 'BP') {
      gene_names_annotated <- AnnotationDbi::select(GOspecies, gene_names, "GO", "SYMBOL")
      gene_names_annotated_DT <- as.data.table(gene_names_annotated)
      Gene <- gene_names_annotated_DT[gene_names_annotated_DT$ONTOLOGY == 'BP']
      gene <- as.data.frame(Gene)
      terms <- AnnotationDbi::select(GO.db, as.character(gene[,2]), "TERM", "GOID")
      symbol_GOID_terms <- cbind(gene$SYMBOL, terms)
      unique(symbol_GOID_terms)
    } else if (deparse(substitute(GOcategory)) == 'CC') {
      gene_names_annotated <- AnnotationDbi::select(GOspecies, gene_names, "GO", "SYMBOL")
      gene_names_annotated_DT <- as.data.table(gene_names_annotated)
      Gene <- gene_names_annotated_DT[gene_names_annotated_DT$ONTOLOGY == 'CC']
      gene <- as.data.frame(Gene)
      terms <- AnnotationDbi::select(GO.db, as.character(gene[,2]), "TERM", "GOID")
      symbol_GOID_terms <- cbind(gene$SYMBOL, terms)
      unique(symbol_GOID_terms)
    } else if (deparse(substitute(GOcategory)) == 'MF') {
      gene_names_annotated <- AnnotationDbi::select(GOspecies, gene_names, "GO", "SYMBOL")
      gene_names_annotated_DT <- as.data.table(gene_names_annotated)
      Gene <- gene_names_annotated_DT[gene_names_annotated_DT$ONTOLOGY == 'MF']
      gene <- as.data.frame(Gene)
      terms <- AnnotationDbi::select(GO.db, as.character(gene[,2]), "TERM", "GOID")
      symbol_GOID_terms <- cbind(gene$SYMBOL, terms)
      unique(symbol_GOID_terms)
    } else
      message("Not a valid GO category.  Must be either BP, CC, or MF.")

  }

}
