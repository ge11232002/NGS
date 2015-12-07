\name{getPairedReads}
\alias{getPairedReads}

\title{
  Get the filenames of paired reads
}
\description{
  Given the filenames of the first reads fastq files, 
  get the filenames of the paired reads fastq if exist.
}
\usage{
  getPairedReads(read1)
}

\arguments{
  \item{read1}{
    The filenames of the first reads.
  }
}
\details{
  The current supported file suffixes are "R1_001.fastq.gz", "R1.fastq.gz",
  "R1.fastq", "F3.csfasta" and "F3.csfasta".
  The corresponding mapping paired reads have suffixes of 
  "R2_001.fastq.gz", "R2.fastq.gz",
  "R2.fastq", "F5-BC.csfasta" and "F5-RNA.csfasta".
  }
\value{
  If the paired reads exist, the filenames are returned.
  Otherwise, NULL is returned.
}

\author{
  Ge Tan
}

\examples{
  x <- c("nanocage_ACAGAT_carp_embryo_R1.fastq.gz",
         "nanocage_CACGAT_carp_hk_control_R1.fastq.gz",
         "nanocage_CACTGA_carp_hk_Tborreli1_R1.fastq.gz",
         "nanocage_CTGACG_carp_hk_Tborreli2_R1.fastq.gz")
  getPairedReads(x)
}
