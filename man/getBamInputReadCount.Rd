\name{getBamInputReadCount}
\alias{getBamInputReadCount}

\title{
  getBamInputReadCount
}
\description{
  Fetch the number of input reads from bam files
}
\usage{
  getBamInputReadCount(bamFiles)
}

\arguments{
  \item{bamFiles}{
    The filenames of input bam files.
  }
}
\details{
  This function reads the header of bam files and searches for the comment line
  "INPUTREADCOUNT", which can be added by \code{\link{addInputReadCounts}}.
}
\value{
  The number of input reads are returned. If the "INPUTREADCOUNT" doesn't exist,
  NA is returned.
}

\author{
  Ge Tan
}

\seealso{
  \code{\link{addInputReadCounts}}
}
\examples{
  \donttest{
  ## This example is not tested because it requires external software "samtools"
  bamfile <- system.file("extdata", "ex1.bam", package="Rsamtools",
                         mustWork=TRUE)
  tempBamFile <- tempfile(fileext=".bam")
  file.copy(bamfile, tempBamFile)
  reads1 <- file.path(system.file("extdata", package="NGS"), "fastq",
                    "nanocage_ACAGAT_carp_embryo_R1.fastq")
  addInputReadCounts(bamFiles=tempBamFile, reads1s=reads1, binary="samtools")
  
  getBamInputReadCount(tempBamFile)
  file.remove(tempBamFile)
  }
}
