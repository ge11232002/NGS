### -----------------------------------------------------------------
### get the fn of the paired reads
### Exported!
getPairedReads <- function(reads1){
  filePairs <- data.frame(reads1=c("R1_001.fastq.gz", "R1.fastq.gz", 
                                   "R1.fastq", "F3.csfasta", "F3.csfasta",
                                   "READ1.fq.gz"),
                          reads2=c("R2_001.fastq.gz", "R2.fastq.gz", 
                                   "R2.fastq", "F5-BC.csfasta", 
                                   "F5-RNA.csfasta",
                                   "READ2.fq.gz"))
  for (i in 1:nrow(filePairs)){
    x <- filePairs$reads1[i]
    if (all(grepPatternList(x, reads1))){
      reads2 <- sub(x, filePairs$reads2[i], reads1)
      if(all(file.exists(reads2))){
        return(reads2)
      }
    }
  }
  return(NULL)
}

### -----------------------------------------------------------------
### count the number of records in fastq files
### Not Exported!
countReadsInFastq <- function(filepath){
  ans <- sapply(lapply(filepath, fastq.geometry), "[", 1)
}

### -----------------------------------------------------------------
### return a list of qualities from shorReadQ class 
### Not Exported!
getQualities <- function(shortReadQ){
  qual <- as(quality(shortReadQ), "matrix")
  qList <- lapply(as.list(as.data.frame(t(qual))), na.omit)
  #qList <- apply(qual, 1, na.omit)
  return(qList)
}

### -----------------------------------------------------------------
### return a list of qualities from shorReadQ class 
### Not Exported!
getLastGoodBasePos <- function(qualities, minTailQuality){
  ## do not search for low qualities in the first 5 bases!!!!
  qualities[1:5] <- minTailQuality 
  lastGood1 <- match(TRUE, runmean(qualities, 4, align="left", alg="fast") <
                       minTailQuality, nomatch=length(qualities)+1) - 1
  return(lastGood1)
}

### -----------------------------------------------------------------
### trimFastq: trim left, right, low quality tail
### Exported!
trimFastqOneRead <- function(reads1, outputDir=".", paired=TRUE, nReads=-1L, 
                             minTailQuality=20, trimLeft=0L, trimRight=0L,
                             mc.cores=getThreads()){
  job <- my.jobStart(paste("start", basename(reads1)))
  
  nTotalReads <- countReadsInFastq(reads1)
  my.writeElapsed(job, status="reads1 counted")
  nTotalReads <- setNames(nTotalReads, basename(reads1))
  if(paired){
    reads2 <- getPairedReads(reads1)
    if(is.null(reads2)){
      stop("The paired reads do not exist!")
    }
    nTotalReads2 <- countReadsInFastq(reads2)
    my.writeElapsed(job, status="reads2 counted")
    nTotalReads2 <- setNames(nTotalReads2, basename(reads2))
    if(any(nTotalReads != nTotalReads2)){
      stop("The number of reads in first reads is differrent 
           from paired reads: ", basename(reads1)[nTotalReads != nTotalReads2])
    }
  }else{
    reads2 <- NULL
  }
  isPaired <- paired
  
  # the numebr of reads can be fewer than  1e6
  nYield <- min(1e6, nTotalReads)
  
  fqs1 <- FastqStreamer(reads1, nYield)
  on.exit(close(fqs1))
  if(isPaired){
    fqs2 <- FastqStreamer(reads2, nYield)
    on.exit(close(fqs2))
  }
  
  ## If we only sample a fixed number of nReads
  if(nReads > 0 && nReads < nTotalReads){ 
    ## nReads is set to -1L by default!
    idx <- round(seq(from=1, to=nYield, 
                     length=(round(nYield/nTotalReads * nReads))))
  }else{
    idx <- NULL
  }
  
  ## Read the fastq files
  while(length(x1 <- yield(fqs1))){
    ## Subset the reads
    if (!is.null(idx)){
      if (length(x1) < nYield){
        ## when reading the last we want to truncate idx
        idx <- idx[idx <= length(x1)]
      }
      x1 <- x1[idx]
    }
    if(isPaired){
      x2 <- yield(fqs2)
      if(!is.null(idx)){
        x2 <- x2[idx]
      }
    }
    endPos1 <- width(x1)
    ## Remove tails
    if(!is.null(minTailQuality)){
      lastGoodBase1 <- unlist(my.mclapply(getQualities(x1), getLastGoodBasePos,
                                          minTailQuality,
                                          mc.cores=mc.cores))
      my.writeElapsed(job, status="getEndPos for reads1")
      use <- lastGoodBase1 >= endPos1/2
      if(isPaired){
        endPos2 <- width(x2)
        lastGoodBase2 <- unlist(my.mclapply(getQualities(x2), 
                                            getLastGoodBasePos,
                                            minTailQuality,
                                            mc.cores=mc.cores))
        my.writeElapsed(job, status="getEndPos for reads2")
        use2 <- lastGoodBase2 >= endPos2/2
        use <- use & use2
        endPos2[use] <- lastGoodBase2[use]
      }
      endPos1[use] <- lastGoodBase1[use]
    }
    if(trimRight > 0){
      endPos1 <- pmin(endPos1, width(x1) - trimRight)
      if(isPaired){
        endPos2 <- pmin(endPos2, width(x2), trimRight)
      }
    }
    target1 <- sub("\\.gz$", "", basename(reads1))
    writeFastq(narrow(x1, start=trimLeft+1, end=endPos1),
               file=file.path(outputDir, target1), mode="a", compress=FALSE)
    if(isPaired){
      target2 <- sub("\\.gz$", "", basename(reads2))
      writeFastq(narrow(x2, start=trimLeft+1, end=endPos2),
                 file=file.path(outputDir, target2), mode="a", compress=FALSE)
    }
  }
  my.writeElapsed(job)
  my.write("file size: ", file.size(reads1)/1e6, " MB")
  if(isPaired){
    my.write("file size: ", file.size(reads2)/1e6, " MB")
  }

  return(target1)
}

trimFastq <- function(reads1, outputDir=".", paired=TRUE, nReads=-1L, 
                      minTailQuality=20, trimLeft=0L, trimRight=0L,
                      mc.cores=getThreads()){
  ans <- sapply(reads1, trimFastqOneRead, outputDir=outputDir,
                paired=paired, nReads=nReads, minTailQuality=minTailQuality,
                trimLeft=trimLeft, trimRight=trimRight,
                mc.cores=mc.cores)
  return(ans)
}

### -----------------------------------------------------------------
### FastUniq: wrapper for De Novo Duplicates Removal Tool for Paired Short Reads
### This only works for paired-end reads and doesn't require reference genome.
### Get rid of duplicates introduced by PCR amplification.
### Exported!
fastUniq <- function(reads1, outputDir=".", binary="fastuniq"){
  if(length(file_ext(reads1)) != 1L){
    stop("The reads1 have different file extensions!")
  }
  
  reads2 <- getPairedReads(reads1)
  if(is.null(reads2)){
    stop("fastUniq only works on paired-end reads!")
  }
  
  if(isGzipped(reads1[1])){
    reads1Temp <- gunzip(filename=reads1, temporary=TRUE)
    reads2Temp <- gunzip(filename=reads2, temporary=TRUE)
  }
  for(i in 1:length(reads1)){
    inputListFile <- tempfile()
    writeLines(c(reads1[i], reads2[i]), con=inputListFile)
    args <- c(paste("-i", inputListFile), 
              paste("-o", file.path(outputDir, sub("\\.gz$", "", 
                                                   basename(reads1[i])))),
              paste("-p", file.path(outputDir, sub("\\.gz$", "",
                                                   basename(reads2[i]))))
    )
    system2(command=binary, args=args)
    unlink(inputListFile)
  }
  return(file.path(outputDir, sub("\\.gz$", "", basename(reads1))))
}

### -----------------------------------------------------------------
### barcodeSplitter: split the fastq file with different barcodes,
### optionally check for must-have sequence at certail position
### Give up. Use TagDust.
### 
barcodeSplitter <- function(reads1, barcodes, mustSeq=NULL, start=NULL,
                            nthread=getThreads()){
  if(!isConstant(nchar(barcodes))){
    stop("The barcodes must have same length!")
  }
  if(xor(is.null(mustSeq), is.null(start))){
    stop("The mustSeq and start must be used or not used together!")
  }
  barcodeWidth <- nchar(barcodes)[1]
  reads2 <- getPairedReads(reads1)
  if(is.null(reads2)){
    message("Single end reads!")
    paired <- FALSE
  }else{
    message("Paired end reads!")
    paired <- TRUE
    fastq2 <- readFastq(reads2)
  }
  fastq1 <- readFastq(reads1)
  
  ## filter out reads without mustSeq
  if(!is.null(mustSeq)){
    mustSubSeqs <- subseq(sread(fastq1), start=start, width=nchar(mustSeq))
    similarityScores <- stringsim(a=mustSeq, b=as.character(mustSubSeqs),
                                  method="jw", nthread=nthread)
    distanceScores <- stringdist(a=mustSeq, b=as.character(mustSubSeqs),
                                  method="osa", nthread=nthread)
  }
  
  ##
  frontSeqs <- subseq(sread(fastq1), start=1, width=barcodeWidth)
  similarityScores <- sapply(barcodes, stringsim, b=as.character(frontSeqs),
                             method="osa", nthread=nthread)
}