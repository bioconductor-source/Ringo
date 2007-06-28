\name{compute.gc}
\alias{compute.gc}
\alias{gccontent}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute the GC content of DNA and probe sequences}
\description{
  Simple auxiliary function to compute the GC content of a given set of
  DNA sequences, such as microarray probe sequences.
}
\usage{
compute.gc(probe.sequences, digits = 2)
}
\arguments{
  \item{probe.sequences}{character vector of DNA or probe sequences of
    which the GC content is to be computed}
  \item{digits}{integer specifying the desired precision}
}
\value{
  a numeric vector with sequence-wise GC contents; the names of this
  vector are the names of the supplied \code{probe.sequences}.
}
\author{Joern Toedling}
\seealso{Function \code{basecontent} in package \code{matchprobes} for a
  more general function to compute base occurrence in sequences}
\examples{
  ex.seqs <- c("gattaca", "GGGNTT", "ggAtT", "tata","gcccg")
  names(ex.seqs) <- paste("sequence",1:5,sep="")
  compute.gc(ex.seqs)
}
\keyword{manip}
