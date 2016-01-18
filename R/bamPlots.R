### -----------------------------------------------------------------
### makeAlignmentCountBarPlot from bam files
### Exported!
makeAlignmentCountBarPlot <- function(bamFiles,
                                      txtFile="read-alignment-statistics.txt",
                                      pairedMode=c("paired", "first", "second"),
                                      mc.cores=getThreads()){
  pairedMode <- match.arg(pairedMode)
  mc.cores <- min(length(bamFiles), mc.cores)
  mmBAM <- mclapply(bamFiles, getBamMultiMatching, pairedMode=pairedMode,
                    mc.cores=mc.cores)
  if(is.null(names(bamFiles))){
    samples <- basename(bamFiles)
  }else{
    samples <- names(bamFiles)
  }
  names(mmBAM) <- samples
  mmValues <- sort(as.integer(unique(unlist(lapply(mmBAM, names)))))
  mmCounts <- matrix(0, nrow=length(samples), ncol=length(mmValues),
                     dimnames=list(samples, mmValues))
  for(sm in samples){
    mm <- mmBAM[[sm]]
    mmCounts[sm, names(mm)] <- mm
  }
  my.write.table(mmCounts, file=txtFile, head="Sample")
  
  ## Plot the results in barplot
  multiCount <- as.integer(colnames(mmCounts))
  isSmall <- multiCount <= 3
  if(any(!isSmall)){
    mmCounts <- cbind(mmCounts[ , isSmall, drop=FALSE],
                      ">3"=rowSums(mmCounts[ , !isSmall, drop=FALSE]))
  }
  multiCountColors <- c("0 hit(s)"="gray", "1 hit(s)"="blue", "2 hit(s)"="cyan",
                        "3 hit(s)"="green", ">3 hit(s)"="orange")
  colnames(mmCounts) <- paste(colnames(mmCounts), "hit(s)")
  stopifnot(colnames(mmCounts) %in% names(multiCountColors))
  par(mar=c(12, 4.1, 4.1, 2.1))
  mmCounts <- t(apply(mmCounts,1,rev))
  barplot(t(mmCounts/1e6), las=2, ylab="Counts [Mio]", 
          main="total alignments", legend.text=TRUE, border=NA,
          col=multiCountColors[colnames(mmCounts)], 
          xlim=c(0, nrow(mmCounts) +5))
}

### -----------------------------------------------------------------
### makeFragLengthHistPlot:  histogram of fragment lengths
### Exported!
makeFragLengthHistPlot <- function(bamFiles, fragSizeMax=1000L,
                                   mc.cores=getThreads()){
  ## First make sure all of bamFiles are paired-end
  paired <- sapply(bamFiles, testPairedEndBam)
  if(!all(paired)){
    stop("The fragment length can only be estimated from paired-end reads: ",
         bamFiles[!paired])
  }
  if(is.null(names(bamFiles))){
    samples <- basename(bamFiles)
  }else{
    samples <- names(bamFiles)
  }
  
  ## load the reads
  mc.cores <- min(mc.cores, length(bamFiles))
  ### The strand is not important here, so we don't specify the strandMode
  reads <- mclapply(bamFiles, readGAlignmentPairs, use.names=TRUE, 
                    mc.cores=mc.cores)
  widths <- list()
  for(i in 1:length(reads)){
    widths[[samples[i]]] <- 
      shrinkToRange(width(unlist(range(grglist(reads[[i]])))), 
                    theRange=c(-0.5, fragSizeMax+0.5))
  }
  toPlot <- data.frame(samples=rep(names(widths), sapply(widths, length)),
                       "fragmentSize"=unlist(widths))
  ggplot(toPlot, aes(x=fragmentSize, fill=samples, colour=samples)) +
    geom_density(alpha = 0.1) + theme_bw() + xlab("Fragment size")
}

### -----------------------------------------------------------------
### makeSegmentCountHistPlots: histogram of aligned segments per read
### Exported!
makeSegmentCountHistPlots <- function(bamFiles,
                                      pairedMode=c("paired", "first",
                                                   "second"),
                                      mc.cores=getThreads()){
  pairedMode <- match.arg(pairedMode)
  mc.cores <- min(length(bamFiles), mc.cores)
  
  paired <- mclapply(bamFiles, testPairedEndBam, mc.cores=mc.cores)
  if(length(unique(unlist(paired))) != 1L){
    stop("All the input bam files have to be both paired-end or single-end")
  }
  
  ## load the bam
  paired <- unique(unlist(paired))
  if(paired){
    param <- ScanBamParam()
    if(pairedMode == "paired"){
      ### The strand is not important here, so we don't specify the strandMode
      reads <- mclapply(bamFiles, readGAlignmentPairs, mc.cores=mc.cores)
    }else if(pairedMode == "first"){
      bamFlag(param) <- scanBamFlag(isFirstMateRead=TRUE, isUnmappedQuery=FALSE)
      reads <- mclapply(bamFiles, readGAlignments,
                        param=param, mc.cores=mc.cores)
    }else if(pairedMode == "second"){
      bamFlag(param) <- scanBamFlag(isSecondMateRead=TRUE,
                                    isUnmappedQuery=FALSE)
      reads <- mclapply(bamFiles, readGAlignments,
                        param=param, mc.cores=mc.cores)
    }
  }else{
    reads <- mclapply(bamFiles, readGAlignments,
                      mc.cores=mc.cores)
  }
  if(is.null(names(bamFiles))){
    samples <- basename(bamFiles)
  }else{
    samples <- names(bamFiles)
  }
  names(reads) <- samples
  
  segmentCounts <- list()
  for(i in 1:length(reads)){
    segmentCounts[[samples[i]]] <- 
      shrinkToRange(njunc(reads[[i]])+1L, theRange=c(0.5, 10.5))
  }
  toPlot <- data.frame(samples=rep(names(segmentCounts),
                                   sapply(segmentCounts, length)),
                       "fragmentSize"=unlist(segmentCounts))
  ggplot(toPlot, aes(x=fragmentSize, fill=samples, colour=samples)) +
    geom_histogram(binwidth=1, position = "dodge") + theme_bw() + 
    xlab("# segments in alignment")
}
