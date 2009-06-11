\name{plot.cher}
\alias{plot.cher}
\alias{plotCher}
\alias{cherPlot}
\alias{plot,cher,missing-method}
\alias{plot,cher,ExpressionSet-method}
\title{Plot identified Chers}
\description{
  Function for plotting identified \emph{Cher}s (ChIP-enriched regions).
}
\usage{
\S4method{plot}{cher,ExpressionSet}(x, y, probeAnno, samples=NULL, extent = 1000, gff = NULL, ...)
}
\arguments{
  \item{x}{object of class \code{cher}}
  \item{y}{data object of class \code{ExpressionSet} that was used for
    function \code{findChersOnSmoothed}}
  \item{probeAnno}{object of class \code{probeAnno}
    holding the reporter/probe to genome mappings}
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
\seealso{\code{\link{chipAlongChrom}}, \code{\link{cher-class}}}
%\examples{}
\keyword{hplot}
