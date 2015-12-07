\name{grepPatternList}
\alias{grepPatternList}

\title{
  Some improved utils functions
}
\description{
  Some base functions are customised for our use.
}
\usage{
  grepPatternList(patternList, x, combine=c("or", "and"), ...)
}

\arguments{
  \item{patternList}{
    The \code{character} string of pattern list to grep.
  }
  \item{x}{
    See the argument \code{x} in \link{\code{grep}} method.
  }
  \item{combine}{
    When it's "or", for each element of x, a TRUE value will be returned
    if it matches to at least one of the pattern list.
    When it's "and", a TRUE value will be returned only if it matches to all of
    the pattern list.
  }
  \item{...}{
    Other arguments passed into \code{grep} method.
  }
}

\value{
  A logical vector of length(x).
}

\author{
  Ge Tan
}

\examples{
  patternList <- c("R1_001.fastq.gz", "R1.fastq.gz", "R1.fastq",
                   "F3.csfasta")
  x <- c("nanocage_ACAGAT_carp_embryo_R1.fastq.gz",
         "nanocage_CACGAT_carp_hk_control_R1.fastq.gz",
         "nanocage_CACTGA_carp_hk_Tborreli1_R1.fastq.gz",
         "nanocage_CTGACG_carp_hk_Tborreli2_R1.fastq.gz")
  grepPatternList(patternList, x, combine="or")
  grepPatternList(patternList, x, combine="and")
}

