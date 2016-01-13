\name{addInputReadCounts}
\alias{addInputReadCounts}

\title{
  addInputReadCounts
}
\description{
  Add the number of input reads as a comment to the header of bamFile.
}
\usage{
  addInputReadCounts(inputs, reads1s, binary = "samtools")
}

\arguments{
  \item{inputs}{
    The filenames of input bam files.
  }
  \item{reads1s}{
    The filenames of the corresponding input fastq files.
  }
  \item{binary}{
    The name/filename of the binary "samtools" to call.
  }
}

\value{
  The filenames of input bam files.
}

\author{
  Ge Tan
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
  addInputReadCounts(inputs=tempBamFile, reads1s=reads1)                  
  file.remove(tempBamFile)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line