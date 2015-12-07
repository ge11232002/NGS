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
trimFastq <- function(reads1, paired=FALSE){
  nTotalReads <- countReadsInFastq(reads1)
  if(paired){
    reads2 <- getPairedReads(reads1)
    if(is.null(reads2)){
      stop("The paired reads do not exist!")
    }
  }
  # the numebr of reads can be fewer than  1e6
  nYield <- min(1e6, nTotalReads)
  
  
}
