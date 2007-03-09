\name{plot.peak}
\alias{plot.peak}
\alias{plotPeak}
\alias{peakPlot}
\title{Plot identified peaks}
\description{
  Function for plotting identified peaks
}
\usage{
\method{plot}{peak}(x, dat, probeAnno, samples=NULL, extent = 1000, gff = NULL, ...)
}
\arguments{
  \item{x}{object of class \code{peak}}
  \item{dat}{data object of class \code{ExpressionSet} that was used for
  function \code{findPeaksOnSmooth}}
  \item{probeAnno}{environment holding probe to genome mappings}
  \item{samples}{which samples to plot, either a numeric vector of
    entries in 1 to \code{ncol(dat)}, or character vector with entries
    in \code{sampleNames(dat)} or \code{NULL} meaning plot the levels
    from all samples in the ExpressionSet}
  \item{extent}{integer; how many base pairs to the left and right
    should the plotted genomic region be extended}
  \item{gff}{data frame with gene/transcript annotation}
  \item{\dots}{further arguments passed on to function
    \code{chipAlongChrom}}
}
\value{
  called for generating the plot; \code{invisible(NULL)}
}
\author{Joern Toedling}
\seealso{\code{\link{chipAlongChrom}}}
\examples{
 \dontrun{
  plot(p, smoothX, probeAnno=exProbeAnno, gff=exGFF)
 }
}
\keyword{hplot}
