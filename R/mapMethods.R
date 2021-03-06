mapBowtie2 <- function(){
  
}

mapBowtie <- function(){
  
}

mapBWA <- function(){
  
}

### -----------------------------------------------------------------
### BWA-MEM algorithm for 70bp-1Mbp query sequences
### Exported!
bwaMem <- function(reads1, output, ref, opt="", mc.cores=getThreads(),
                   binaryBWA="bwa", binarySamtools="samtools"){
  if(length(reads1) != 1L || length(output) != 1L){
    stop("Please provide one sample each time!")
  }
  if(file_ext(output) != "bam"){
    stop("Please specify the output name with bam as file extension.")
  }
  reads2 <- getPairedReads(reads1)
  if(is.null(reads2)){
    paired <- FALSE
    reads2 <- ""
  }else{
    paired <- TRUE
  }
  message("Paired? ", paired)
  args <- c("mem", "-t", mc.cores, opt, ref, reads1, reads2,
            "|", binarySamtools, "view -b - >", output)
  ## run bwa
  res <- system2(command=binaryBWA, args=args)
  stopifnot(res == 0L)
  ### res: 0 success
  
  ## add input read count
  message("Add input read counts.")
  addInputReadCounts(output, reads1)
  
  ## sort and index bam
  message("Sort and index bam files.")
  tempBam <- tempfile(fileext=".bam")
  tempBam <- sortBam(output, sub("\\.bam$", "", tempBam, ignore.case=TRUE))
  
  stopifnot(file.copy(tempBam, output, overwrite=TRUE))
  stopifnot(file.remove(tempBam))
  indexBam(output)
  
  ## bam to bigwig files
  message("bam to bigwig conversion.")
  bam2bigwig(output)
  
  invisible(output)
}

mapTophat <- function(){
  
}
