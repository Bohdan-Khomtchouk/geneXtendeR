geneXtender_p <- function(upstream, organism) {
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
  levels(geneXtender.file$seqid) <- gsub("chr", "", levels(geneXtender.file$seqid))
  levels(geneXtender.file$seqid) <- gsub("X", "100", levels(geneXtender.file$seqid))
  levels(geneXtender.file$seqid) <- gsub("Y", "200", levels(geneXtender.file$seqid))
  levels(geneXtender.file$seqid) <- gsub("MT", "300", levels(geneXtender.file$seqid))
  levels(geneXtender.file$seqid) <- gsub("MtDNA", "300", levels(geneXtender.file$seqid))
  levels(geneXtender.file$seqid) <- gsub("M", "300", levels(geneXtender.file$seqid))
  levels(geneXtender.file$seqid) <- gsub("Mito", "300", levels(geneXtender.file$seqid))
  levels(geneXtender.file$seqid) <- gsub("I", "1", levels(geneXtender.file$seqid))
  levels(geneXtender.file$seqid) <- gsub("II", "2", levels(geneXtender.file$seqid))
  levels(geneXtender.file$seqid) <- gsub("III", "3", levels(geneXtender.file$seqid))
  levels(geneXtender.file$seqid) <- gsub("IV", "4", levels(geneXtender.file$seqid))
  levels(geneXtender.file$seqid) <- gsub("V", "5", levels(geneXtender.file$seqid))
  levels(geneXtender.file$seqid) <- gsub("VI", "6", levels(geneXtender.file$seqid))
  levels(geneXtender.file$seqid) <- gsub("VII", "7", levels(geneXtender.file$seqid))
  levels(geneXtender.file$seqid) <- gsub("VIII", "8", levels(geneXtender.file$seqid))
  levels(geneXtender.file$seqid) <- gsub("IX", "9", levels(geneXtender.file$seqid))
  levels(geneXtender.file$seqid) <- gsub("XI", "11", levels(geneXtender.file$seqid))
  levels(geneXtender.file$seqid) <- gsub("XII", "12", levels(geneXtender.file$seqid))
  levels(geneXtender.file$seqid) <- gsub("XIII", "13", levels(geneXtender.file$seqid))
  levels(geneXtender.file$seqid) <- gsub("XIV", "14", levels(geneXtender.file$seqid))
  levels(geneXtender.file$seqid) <- gsub("XV", "15", levels(geneXtender.file$seqid))
  levels(geneXtender.file$seqid) <- gsub("XVI", "16", levels(geneXtender.file$seqid))
  geneXtender.file$seqid = as.numeric(as.character(geneXtender.file$seqid))
  geneXtender.sorted <- dplyr::arrange(geneXtender.file, as.numeric(seqid), start)
  write.table(geneXtender.sorted, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE, sprintf("geneXtender_gtf_%s.bed", upstream))
  return(data.table::as.data.table(geneXtender.sorted))
}


annotate_n <- function(organism, extension, n=2, aggregation=TRUE) {
  if(!file.exists("peaks.txt")){
    message("Please run peaksInput() function first!  See ?peaksInput for more information")
  } else if(n < 1) {
    message("Please run a valid distance of gene's away with n > 0")
  }else {
    oopts = options(warn=-1)
    on.exit(options(oopts))
    
    genes <- geneXtender_p(extension, organism)
    
    run3 <- function(f1, f2, peakslist) {
      .C("annotate", f1, f2, peakslist)[[3]]
    }
    
    linelen <- "                                                                                                    "
    peaksArray<-rep(linelen,500000)
    sapply(extension, geneXtender_p, organism=organism)
    onegxFile <- sprintf("geneXtender_gtf_%s.bed", extension)
    onecmd2 <- run3(f1 = "peaks.txt", f2 = onegxFile, peakslist = peaksArray) 
    onecmd3 <- onecmd2[onecmd2 != linelen]
    
    rdt <- data.table::transpose(data.table::as.data.table(strsplit(onecmd3, "\t"))) #Formatting
    colnames(rdt) <- c("Chromosome", "Peak-Start", "Peak-End", "rm", "Gene-Start", "Gene-End", "Gene-ID", "Gene-Name", "Distance-of-Gene-to-Nearest-Peak")
    rdt <- rdt[,rm:= NULL]
    for (j in c("Chromosome", "Peak-Start", "Peak-End", "Gene-Start", "Gene-End", "Distance-of-Gene-to-Nearest-Peak")) data.table::set(rdt, j=j, value = as.numeric(rdt[[j]]))
    
    indices <- sapply(rdt$`Gene-ID`, function(x) which(x == genes$gene_id))
    
    rdj <- rdt[,c("seqid", "start", "end", "gene_id", "gene_name") := genes[indices+(n-1)]]
    
    err <- which(rdj$Chromosome != rdj$seqid)
    if (length(err) > 0) {
      message(sprintf("Warning! There were %s peaks that could not find the next gene on their respective chromosomes %s genes away", errlen, n))
      for(col in c("seqid", "start", "end", "gene_id", "gene_name")) set(rdj, i=err, j=col, value=NA) #cleanup
    }
    
    setnames(rdj, c("seqid", "start", "end", "gene_id", "gene_name"), c("Chromosome-for-Gene-N-away", "Start-Gene-N-away", "End-Gene-N-away", "Gene-ID-N-away", "Gene-Name-N-away"))
    return(rdt)
  }
}