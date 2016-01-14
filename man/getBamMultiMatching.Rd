\name{getBamMultiMatching}
\alias{getBamMultiMatching}

\title{
  getBamMultiMatching
}
\description{
  Get the numbers of multi-hits for bam file
}
\usage{
  getBamMultiMatching(bamFile, pairedMode = c("paired", "first", "second"))
}

\arguments{
  \item{bamFile}{
    The filenames of input bam file.
  }
  \item{pairedMode}{
    The \code{pairedMode} can be "paired", "first", "second".
    "paired" will load proper pairs; "first" will load all the first reads;
    "second" will load all mapped mate reads.
    This argument is ignored when it's single-end read.
  }
}

\value{
  It returns a integer vector that holds the numbers of unmapped reads 
  (If the input read count cannot be retrieved from bam file), 
  the reads that have exactly one hit, the reads that have exactly two hits,
  and so on.
}

\author{
  Ge Tan
}

\seealso{
  \code{\link{getBamInputReadCount}}
}
\examples{
  \dontrun{
  ## This example is not tested because it requires external software "samtools"
  bamfile <- system.file("extdata", "ex1.bam", package="Rsamtools",
                         mustWork=TRUE)
  tempBamFile <- tempfile(fileext=".bam")
  file.copy(bamfile, tempBamFile)
  reads1 <- file.path(system.file("extdata", package="NGS"), "fastq",
                    "nanocage_ACAGAT_carp_embryo_R1.fastq")
  addInputReadCounts(bamFiles=tempBamFile, reads1s=reads1, binary="samtools")
  
  getBamMultiMatching(tempBamFile)
  file.remove(tempBamFile)
  }
}
