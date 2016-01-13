\name{fastUniq}
\alias{fastUniq}

\title{
  Duplicated reads removal
}
\description{
  Wrapper function of FastUniq: 
  A Fast De Novo Duplicates Removal Tool for Paired Short Reads.
}
\usage{
  fastUniq(reads1, outputDir=".", binary="fastuniq")
}

\arguments{
  \item{reads1}{
    The filenames of first reads.
  }
  \item{outputDir}{
    Where to put the output.
  }
  \item{binary}{
    The name/filename of the binary to call.
  }
}
\details{
  The original FastUniq works for both Fastq and Fasta files.
  We only use it for fastq files.
  
  Other possible tools for removing duplicates from fastq files:
  1. fastx_collapser in the FASTX-Toolkit; 2. rmdup in the SAMtools package;
  3. MarkDuplicates in the Picard toolkit.
}
\value{
  The filenames of generated first reads.
}
\references{
Xu, H., Luo, X., Qian, J., Pang, X., Song, J., Qian, G., Chen, J., and Chen, S. (2012). FastUniq: A Fast De Novo Duplicates Removal Tool for Paired Short Reads. PLoS One 7.
}
\author{
  Ge Tan
}

\note{
  Usually it's not recommended to remove duplicates in RNA-Seq dataset.
}

\examples{
  \donttest{
  ## This example is not tested because it requires external software "FastUniq"
    reads1 <- file.path(system.file("extdata", package="NGS"), "fastq",
                        "nanocage_ACAGAT_carp_embryo_R1.fastq")
    fastUniq(reads1, outputDir=tempdir())
  }
}
