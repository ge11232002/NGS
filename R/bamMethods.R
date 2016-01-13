### -----------------------------------------------------------------
### sortIndexBam: sort and index a bam file
### Exported!!
sortIndexBam <- function(input, output,
                         maxMem="4096M", removeBam=TRUE,
                         mc.cores=getThreads(), binary="samtools"){
  if(length(input) != length(output)){
    stop("The number of input bam names differ from output bam names!")
  }
  for(i in 1:length(input)){
    args <- c("sort",
              paste("-m", maxMem),
              paste("-@", mc.cores),
              paste("-o", output[i]),
              input[i])
    
    system2(command=binary, args=args)
    if(removeBam){
      file.remove(input[i])
    }
    args <- c("index", output[i])
    system2(command=binary, args=args)
  }
  invisible(output)
}

