\name{bwaMem}
\alias{bwaMem}

\title{
  bwaMem
}
\description{
  Wrapper function for bwa aligner with MEM algorithm.
}
\usage{
  bwaMem(reads1, output, ref, opt = "", mc.cores = getThreads(), 
         binaryBWA = "bwa", binarySamtools = "samtools")
}

\arguments{
  \item{reads1}{
    The filename of first reads.
  }
  \item{output}{
    The filename of output bam files.
  }
  \item{ref}{
    The filename of bwa index prefix.
  }
  \item{opt}{
    Additional arguments passed to bwa-mem.
  }
  \item{mc.cores}{
    The number of threads to use. By default,
    we use a internal function  \code{NGS:::getThreads} to decide.
  }
  \item{binaryBWA}{
    The name/filename of the binary "bwa" to call.
  }
  \item{binarySamtools}{
    The name/filename of the binary "samtools" to call.
  }
}
\details{
  The BWA-MEM algorithm is suitable for aligning 70bp-1Mbp query sequences.
  The bam file will be sorted by coordinates.
  A idnex file will be created.
  A corresponding bigwig file will also be created.
}
\value{
  The filename of output bam file. 
}
\references{
  http://bio-bwa.sourceforge.net/bwa.shtml
}
\author{
  Ge Tan
}


\examples{
  \dontrun{
  ## This example is not tested because it requires external software "bwa" and
  ## reads/reference
    bwaMem(reads1="../../fastq/trim_galore/nanocage_ACAGAT_carp_embryo_R1_val_1.fq.gz", 
           output="nanocage_ACAGAT_carp_embryo.bam", 
           ref="/mnt/biggles/data/reference/Cyprinus_carpio/Leiden/cypCar1/Sequence/BWAIndex/genome.fa", 
           opt = "", mc.cores=4)
  }
}

