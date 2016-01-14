\name{trimFastq}
\alias{trimFastq}

\title{
  Trim the fastq files
}
\description{
  Given the argument of trimLeft, trimRight and minTailQuality,
  this function trims the fastq files and generate the new files.
}
\usage{
  trimFastq(reads1, outputDir=".", paired=TRUE, nReads=-1L, 
            minTailQuality=20, trimLeft=0L, trimRight=0L,
            mc.cores=getThreads())
}

\arguments{
  \item{reads1}{
    The \code{character} vector. The filenames of first reads.
  }
  \item{outputDir}{
    When to save the trimmed fastq files.
  }
  \item{paired}{
    Whether paired-end reads or single-end reads.
  }
  \item{nReads}{
    Only sample a certain number \code{nReads} of reads.
    By default, all reads will be processed.
  }
  \item{minTailQuality}{
    The minimal tail quaility to decided where to cut.
    We use a runnning windows of size 4 to calculate the quality at each base,
    and cut off the tail at the first base where the quality is below
    \code{minTailQuality}. By default, it is 20.
    If the length of the trimmed read is below half of original read,
    this read is discarded.
  }
  \item{trimLeft}{
    The number of bases to cut on the left of the read.
  }
  \item{trimRight}{
    The number of bases to cut on the right (tail) of the read.
  }
  \item{mc.cores}{
    The number of threads to use. By default,
    we use a internal function  \code{NGS:::getThreads} to decide.
  }
}

\value{
  The filenames of trimmed first reads are returned.
}

\author{
  Ge Tan
}


\examples{
  reads1 <- file.path(system.file("extdata", package="NGS"), "fastq", 
                      "nanocage_ACAGAT_carp_embryo_R1.fastq")
  trimmedFastqs <- trimFastq(reads1, outputDir=".", paired=TRUE, 
                             nReads=-1L, 
                             minTailQuality=20, 
                             trimLeft=0L, trimRight=0L,
                             mc.cores=1L)
  file.remove(trimmedFastqs, NGS:::getPairedReads(trimmedFastqs))
}

