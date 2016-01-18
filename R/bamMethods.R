### -----------------------------------------------------------------
### sortIndexBam: sort and index a bam file
### Exported!!
sortIndexBam <- function(inputs, outputs,
                         maxMem="4096M", removeBam=TRUE,
                         mc.cores=getThreads(), binary="samtools"){
  if(length(inputs) != length(outputs)){
    stop("The number of input bam names differ from output bam names!")
  }
  for(i in 1:length(inputs)){
    args <- c("sort",
              paste("-m", maxMem),
              paste("-@", mc.cores),
              paste("-o", outputs[i]),
              inputs[i])
    
    system2(command=binary, args=args)
    if(removeBam){
      file.remove(inputs[i])
    }
    args <- c("index", outputs[i])
    system2(command=binary, args=args)
  }
  invisible(outputs)
}

### -----------------------------------------------------------------
### adds the number of reads in reads1 as a comment to the header of bamFile
### Exported!
addInputReadCounts <- function(bamFiles, reads1s, binary="samtools"){
  if(length(bamFiles) != length(reads1s)){
    stop("The number of input bam names differ from input fastq names!")
  }
  for(i in 1:length(bamFiles)){
    bamFile <- bamFiles[i]
    reads1 <- reads1s[i]
    
    # count the reads in the reads1 file
    nReads <- countReadsInFastq(reads1)
    
    # use sammtools view -H bamFile > header.txt
    headerFile <- tempfile(fileext=".header")
    args <- c("view -H", bamFile, ">", headerFile)
    system2(command=binary, args=args)
    
    # append header.txt with a line CO <TAB>INPUTREADCOUNT:NNNN
    cat("@CO\tINPUTREADCOUNT:", format(nReads, scientific=FALSE), "\n", 
        sep="", file=headerFile, append=TRUE)
    
    # use "samtools reheader header.txt bamFile > out.bam" to change the header
    tempBam <- tempfile(fileext=".bam")
    args <- c("reheader", headerFile, bamFile, ">", tempBam)
    system2(command=binary, args=args)
    file.remove(headerFile)
    
    # Move the tempBam back
    file.copy(tempBam, bamFile, overwrite=TRUE)
    file.remove(tempBam)
    
    # recreate the index file if exists
    baiFile <- paste0(bamFile, ".bai")
    if(file.exists(baiFile)){
      indexBam(bamFile)
    }
  }
  invisible(bamFiles)
}

### -----------------------------------------------------------------
### bam2bigwig: convert a bam file into bw file with coverage along the genome
### We load all the mapped reads: no matter it's paired/proper paired or not
### Exported!
bam2bigwig <- function(bamFiles, 
                       bigwigFiles=sub("\\.bam$", ".bw", bamFiles, 
                                       ignore.case=TRUE)
                       ){
  if(length(bamFiles) != length(bigwigFiles)){
    stop("The number of input bam names differ from output bigwig names!")
  }
  for(i in 1:length(bamFiles)){
    job <- my.jobStart(paste("start", basename(bamFiles)))
    gal1 <- readGAlignments(bamFiles[i])
    cov1 <- coverage(gal1)
    export.bw(cov1, bigwigFiles[i])
    my.writeElapsed(job, status="finished writing bigwig")
    rm(gal1, cov1)
  }
  invisible(bigwigFiles)
}

### -----------------------------------------------------------------
### getBamMultiMatching: get the numbers of multi-hits for bam file
### The input bam files have to be sorted and indexed.
### Exported!
getBamMultiMatching <- function(bamFile,
                                pairedMode=c("paired", "first", "second")){
  pairedMode <- match.arg(pairedMode)
  
  baiFile <- paste0(bamFile, ".bai")
  if(!file.exists(baiFile)){
    stop("The input bam file has to be sorted and indexed!")
  }
  
  job <- my.jobStart(paste("bam multimatch", basename(bamFile)))
  param <- ScanBamParam(what="qname")
  bamFlag(param) <- scanBamFlag(isUnmappedQuery=FALSE)
  
  ## test paired-end or not
  paired <- testPairedEndBam(bamFile)
  if(paired){
    bamFlag(param) <- switch(pairedMode,
                             paired=scanBamFlag(isFirstMateRead=TRUE,
                                                isProperPair=TRUE,
                                                isUnmappedQuery=FALSE),
                             first=scanBamFlag(isFirstMateRead=TRUE,
                                               isUnmappedQuery=FALSE),
                             second=scanBamFlag(isSecondMateRead=TRUE,
                                                isUnmappedQuery=FALSE))
  }
  bamReads <- scanBam(bamFile, param=param)[[1]]$qname
  result <- table(bamReads)
  result2 <- table(result)
  
  nReads <- getBamInputReadCount(bamFile)
  if(!is.na(nReads)){
    nReadsUnmapped <- nReads - sum(result)
    result2 <- c("0"=nReadsUnmapped, result2)
    my.writeElapsed(job, "nreads from header loaded")
  }
  my.writeElapsed(job, "done")
  return(result2)
}

### -----------------------------------------------------------------
### getBamInputReadCount: get the number of input reads
### Exported!
getBamInputReadCount <- function(bamFiles){
  ans <- c()
  for(i in 1:length(bamFiles)){
    bamFile <- bamFiles[i]
    header <- scanBamHeader(bamFile)[[1]]
    ## test whether there is any comment line
    if(sum(names(header$text) == "@CO") == 0){
      ans <- c(ans, NA)
    }
    ## get the comment lines
    commentLinesIndex <- which(names(header$text) == "@CO")
    for(i in 1:length(commentLinesIndex)){
      commentLine <- header$text[[commentLinesIndex[i]]]
      searchResults <- grepl("^INPUTREADCOUNT:", commentLine)
      numOfCount <- sum(searchResults)
      if(numOfCount == 0){
        next
      }else if(numOfCount == 1){
        nReads <- as.integer(sub("^INPUTREADCOUNT:", "", commentLine[which(searchResults)]))
        ans <- c(ans, nReads)
      }else{
        stop("more than one INPUTREADCOUNT found in ", bamFile)
      }
    }
  }
  return(ans)
}


