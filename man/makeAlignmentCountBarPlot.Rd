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
  bamFile <- system.file("extdata", "bam", "nanocage_ACAGAT_carp_embryo.bam",
                         package="NGS")
  ### Plot the barplot
  makeAlignmentCountBarPlot(bamFile)
  file.remove("read-alignment-statistics.txt")
  }
}
