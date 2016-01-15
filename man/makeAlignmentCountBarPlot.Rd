\name{makeAlignmentCountBarPlot}
\alias{makeAlignmentCountBarPlot}

\title{
  makeAlignmentCountBarPlot
}
\description{
  Make the plot of Multi-Matching Reported in Bam Files.
}
\usage{
  makeAlignmentCountBarPlot(bamFiles, txtFile="read-alignment-statistics.txt",
                            pairedMode=c("paired", "first", "second"),
                            mc.cores=getThreads())
}

\arguments{
  \item{bamFiles}{
    The filenames of input bam files.
  }
  \item{txtFile}{
    The filename of output for the summary table.
  }
  \item{pairedMode}{
    The \code{pairedMode} can be "paired", "first", "second".
    "paired" will load proper pairs; "first" will load all the first reads;
    "second" will load all mapped mate reads.
    This argument is ignored when it's single-end read.
  }
  \item{mc.cores}{
    The number of threads to use. By default, 
    we use a internal function  \code{NGS:::getThreads} to decide.
  }
}
\details{
  The barplot gives the number of unmapped reads 
  (If the input read count cannot be retrieved from bam file),
  the reads that have exactly one hit, the reads that have exactly two hits,
  and so on.
}
\value{
  The object returned by barplot.
}

\author{
  Ge Tan
}

\seealso{
  \code{\link{getBamMultiMatching}}
}
\examples{
  \donttest{
  ## This example is not tested because it requires external software "samtools"
  ### Prepare the bam files
  bamfile <- system.file("extdata", "ex1.bam", package="Rsamtools",
                         mustWork=TRUE)
  tempBamFile <- tempfile(fileext=".bam")
  file.copy(bamfile, tempBamFile)
  reads1 <- file.path(system.file("extdata", package="NGS"), "fastq",
                    "nanocage_ACAGAT_carp_embryo_R1.fastq")
  addInputReadCounts(bamFiles=tempBamFile, reads1s=reads1, binary="samtools")
  tempBamFile2 <- tempfile(fileext=".bam")
  tempBamFile2 <- sortBam(tempBamFile, 
                     sub("\\\.bam$", "", tempBamFile2, ignore.case=TRUE))
  file.copy(tempBamFile2, tempBamFile, overwrite=TRUE)
  file.remove(tempBamFile2)
  indexBam(tempBamFile)
  
  ### Plot the barplot
  makeAlignmentCountBarPlot(tempBamFile)
  file.remove(read-alignment-statistics.txt)
  file.remove(tempBamFile)
  }
}
