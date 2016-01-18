\name{makeSegmentCountHistPlots}
\alias{makeSegmentCountHistPlots}

\title{
  makeSegmentCountHistPlots
}
\description{
  Plot the histogram of the number of alignments per read.
}
\usage{
  makeSegmentCountHistPlots(bamFiles,
                            pairedMode = c("paired", "first", "second"),
                            mc.cores = getThreads())
}

\arguments{
  \item{bamFiles}{
    The \code{character}() filenames of input bam files.
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
  When it's paired-end reads, the number of alignments is calculated for
  the whole pair (fragment).
}
\value{
  A \code{ggplot} object.
}

\author{
  Ge Tan
}

\seealso{
  \code{\link{makeFragLengthHistPlot}}
}
\examples{
  bamfile <- system.file("extdata", "ex1.bam", package="Rsamtools",
                         mustWork=TRUE)
  makeSegmentCountHistPlots(bamfile, pairedMode="paired", mc.cores=1L)
}
