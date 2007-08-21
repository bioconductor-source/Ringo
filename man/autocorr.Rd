\name{autocor}
\alias{autocorr}
\alias{autocor}
\alias{autocorrelation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function to compute auto-correlation of probe intensities}
\description{
  Function to compute auto-correlation of probe intensities at specified
  offsets from the original positions. 
}
\usage{
autocor(x, probeAnno, chrom = "1", samples = NULL, lag.max = 2000,
lag.step=100, cor.method = "pearson", 
channel = c("red","green","logratio"), verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object either of class \code{ExpressionSet} containing the
    normalized probe intensities or of class \code{RGList} containing
    the raw intensities.}
  \item{probeAnno}{Environment holding chromosomal start positions and
    indices of probes in data matrix.}
  \item{chrom}{character; chromosome to compute the autocorrelation for}
  \item{samples}{which samples of the data to use; if more than 1 for
    each probe the mean intensity over these samples is taken.}
  \item{lag.max}{integer; maximal offset from the original position, the
    auto-correlation is to be computed for.}
  \item{lag.step}{integer; step size of lags between 0 and maximal lag.}
   \item{cor.method}{character; which type of correlation to compute,
    translates to argument \code{method} of function \code{cor}}
  \item{channel}{character; in case \code{x} is an \code{RGList}, which
    channel to plot, either \code{red}, \code{green} or the
    \code{logratio} \code{log2(red)-log2(green)}}
  \item{verbose}{logical; extended output to STDOUT}
}
\details{
  the lags, i.e. the offsets from the original position, the
  auto-correlation is to be computed for, are constructed from the
  given arguments as \code{seq(0, lag.max, by=lag.step)}.
}
\value{
  Object of class \code{autocor.result}, a
  numeric vector of auto-correlation values at the offsets specified in
  argument \code{lags}. The lag values are stored as the \code{names} of
  the vector. Argument \code{chrom} is stored as attribute
  \code{chromosome} of the result.
}
\author{Joern Toedling \email{toedling@ebi.ac.uk}}
\seealso{\code{\link{cor}},\code{\link{plot.autocor.result}}}
\examples{
  exDir <- system.file("exData",package="Ringo")
  load(file.path(exDir,"exampleProbeAnno.rda"))
  load(file.path(exDir,"exampleX.rda"))
  exAc <- autocor(exampleX, probeAnno=exProbeAnno,
                  chrom="9", lag.max=1000)
  plot(exAc)
}
\keyword{manip}
