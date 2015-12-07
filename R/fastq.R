### -----------------------------------------------------------------
### get the fn of the paired reads
### Exported!
getPairedReads <- function(reads1){
  filePairs <- data.frame(reads1=c("R1_001.fastq.gz", "R1.fastq.gz", 
                                   "R1.fastq", "F3.csfasta", "F3.csfasta"),
                          reads2=c("R2_001.fastq.gz", "R2.fastq.gz", 
                                   "R2.fastq", "F5-BC.csfasta", 
                                   "F5-RNA.csfasta"))
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
### trimFastq: trim left, right, low quality tail
### Exported!
trimFastq <- function(reads1, paired=TRUE, nReads=-1L, minTailQuality=20){
  nTotalReads <- countReadsInFastq(reads1)
  if(paired){
    reads2 <- getPairedReads(reads1)
    if(is.null(reads2)){
      stop("The paired reads do not exist!")
    }
  }else{
    reads2 <- NULL
  }
  isPaired <- paired 
  
  # the numebr of reads can be fewer than  1e6
  nYield <- min(1e6, nTotalReads)
  
  fqs1 <- FastqStreamer(reads1, nYield)
  if(isPaired){
    fqs2 <- FastqStreamer(reads2, nYield)
  }
  
  ## If we only sample a fixed number of nReads
  if(config$nReads > 0 && config$nReads < nTotalReads){ 
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
    if (isPaired){
      x2 <- yield(fqs2)
      if (!is.null(idx)){
        x2 <- x2[idx]
      }
    }
    
    ## Remove tails
    if (!is.null(minTailQuality)){
      
    }
  }
}
