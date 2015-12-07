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
    if (my.grep(x, reads1)){
      reads2 <- sub(x, filePairs$reads2[i], reads1)
      if(file.exists(reads2)){
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
