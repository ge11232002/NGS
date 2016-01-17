\name{makeFragLengthHistPlot}
\alias{makeFragLengthHistPlot}

\title{
  makeFragLengthHistPlot
}
\description{
  Make the histogram distribution of fragment lengths estimated from paired-end
  reads.
}
\usage{
  makeFragLengthHistPlot(bamFiles, fragSizeMax = 1000L,
                         mc.cores = getThreads())
}

\arguments{
  \item{bamFiles}{
    The \code{character}() filenames of input bam files.
  }
  \item{fragSizeMax}{
    The \code{integer}(1) of maximal fragment size to plot.
  }
  \item{mc.cores}{
    The number of threads to use. By default,
    we use a internal function  \code{NGS:::getThreads} to decide.
  }
}

\value{
  A \code{ggplot} object.
}

\author{
  Ge Tan
}

\seealso{
  \code{\link{makeAlignmentCountBarPlot}}
}
\examples{
  bamfile <- system.file("extdata", "ex1.bam", package="Rsamtools",
                         mustWork=TRUE)
  makeFragLengthHistPlot(bamfile, fragSizeMax=500L, mc.cores=1L)
}
