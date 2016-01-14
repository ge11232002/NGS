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
addInputReadCounts <- function(inputs, reads1s, binary="samtools"){
  if(length(inputs) != length(reads1s)){
    stop("The number of input bam names differ from input fastq names!")
  }
  for(i in 1:length(inputs)){
    bamFile <- inputs[i]
    reads1 <- reads1s[i]
    
    # count the reads in the reads1 file
    nReads <- countReadsInFastq(reads1)
    
    # use sammtools view -H bamFile > header.txt
    headerFile <- tempfile(fileext="header")
    args <- c("view -H", bamFile, ">", headerFile)
    system2(command=binary, args=args)
    
    # append header.txt with a line CO <TAB>INPUTREADCOUNT:NNNN
    cat("@CO\tINPUTREADCOUNT:", format(nReads, scientific=FALSE), "\n", 
        sep="", file=headerFile, append=TRUE)
    
    # use "samtools reheader header.txt bamFile > out.bam" to change the header
    tempBam <- tempfile(fileext="bam")
    args <- c("reheader", headerFile, bamFile, ">", tempBam)
    system2(command=binary, args=args)
    
    # Move the tempBam back
    file.rename(tempBam, bamFile)
    
    # recreate the index file
    indexBam(bamFile)
    
    file.remove(headerFile)
  }
  invisible(inputs)
}

### -----------------------------------------------------------------
### bam2bigwig: refer to http://rseqc.sourceforge.net/#bam2wig-py
### Exported!
bam2bigwig <- function(){
  
}