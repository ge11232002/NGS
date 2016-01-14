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
    
    # Move the tempBam back
    file.copy(tempBam, bamFile, overwrite=TRUE)
    file.remove(tempBam)
    
    # recreate the index file
    indexBam(bamFile)
    
    file.remove(headerFile)
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

