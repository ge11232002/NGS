\name{bam2bigwig}
\alias{bam2bigwig}

\title{
  bam2bigwig
}
\description{
  Convert a bam file into a bigwig file with coverage.
}
\usage{
  bam2bigwig(bamFiles, bigwigFiles = sub("\\\.bam$", ".bw", bamFiles, 
                                         ignore.case = TRUE))
}

\arguments{
  \item{bamFiles}{
    The \code{character}() filenames of input bam files.
  }
  \item{bigwigFiles}{
    The filenames of output bigwig files.
  }
}
\details{
  All the mapped reads are loaded, no matter they are paired or not,
  properly paired or not.
}
\value{
  The filenames of output bigwig files.
}

\author{
  Ge Tan
}

\examples{
  bamfile <- system.file("extdata", "ex1.bam", package="Rsamtools",
                         mustWork=TRUE)
  bwFile <- tempfile(fileext=".bw")
  bam2bigwig(bamfile, bwFile)
  file.remove(bwFile)
}
