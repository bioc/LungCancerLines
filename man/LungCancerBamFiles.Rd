\name{LungCancerBamFiles}
\alias{LungCancerBamFiles}
\title{
  Get the BAM paths
}
\description{
  Gets a \code{BamFileList} pointing to BAM files containing read
  alignments for the H1993 and H2073 RNA-seq samples. The files are
  the \dQuote{analyzed} BAM files as output by the \code{HTSeqGenie}
  package. 
}
\usage{
LungCancerBamFiles()
}
\details{
  The reads were aligned to genome \code{\link[gmapR]{TP53Genome}},
  using the following parameters:
  \itemize{    
    \item{\code{splicing}: knownGene}
    \item{\code{novelsplicing}: 1}
    \item{\code{indel_penalty}: 1}
    \item{\code{distant_splice_penalty}: 1}
    \item{\code{suboptimal_levels}: 2}
    \item{\code{npaths}: 10}
  }
  Note that the BAM files contain only unique alignments.
}
\value{
  A \code{BamFileList} pointing to two BAM files, one for H1993, one for
  H2073.
}
\author{
  Michael Lawrence
}

\examples{
LungCancerBamFiles()
}
