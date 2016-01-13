\name{sortIndexBam}
\alias{sortIndexBam}

\title{
  sortIndexBam
}
\description{
  sortIndexBam sort the bam file by coordinates and index it.
}
\usage{
  sortIndexBam <- function(inputs, outputs, maxMem="4096M", removeBam=TRUE,
                           mc.cores=getThreads(), binary="samtools")
}

\arguments{
  \item{inputs}{
    The filenames of input bam files.
  }
  \item{outputs}{
    The filenames of output bam files.
  }
  \item{maxMem}{
    The maximal memory per thread to use.
  }
  \item{removeBam}{
    Remove the input bam files or not after sorting and idnexing.
  }
  \item{mc.cores}{
    The number of threads to use. By default, 
    we use a internal function  \code{NGS:::getThreads} to decide.
  }
  \item{binary}{
    The name/filename of the binary "samtools" to call.
  }
}

\value{
  The output is returned.
}

\author{
  Ge Tan
}

\seealso{
  indexBam and sortBam, asBam from Rsamtools
}
\examples{
  \donttest{
  ## This example is not tested because it requires external software "samtools"
  library(Rsamtools)
  bamfile <- system.file("extdata", "ex1.bam", package="Rsamtools",
                           mustWork=TRUE)
  sortIndexBam(inputs=bamfile, outputs="ex1_sorted.bam", removeBam=FALSE)
  }
}