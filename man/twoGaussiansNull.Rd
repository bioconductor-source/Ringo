\name{twoGaussiansNull}
\alias{twoGaussiansNull}
\title{
  Estimate a threshold from Gaussian mixture distribution
}
\description{
  Function to estimate a threshold from Gaussian mixture distribution.
  The data is assumed to follow a mixture of two Gaussian
  distributions. The one Gaussian with the lower mean value is assumed
  to be the null distribution and probe levels are assigned p-values
  based on this null distribution. The threshold is then the minimal
  data value with an adjusted p-value smaller than a specified
  threshold.
}
\usage{
twoGaussiansNull(x, p.adj.method = "BY", max.adj.p = 0.1, var.equal = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{numeric vector of data values}
  \item{p.adj.method}{method for adjusting the p-values for multiple
    testing; must be one of \code{p.adjust.methods}}
  \item{max.adj.p}{which adjusted p-value to use as upper limit for
    estimating the threshold}
  \item{var.equal}{logical; is the variance of the two Gaussians assumed
    to be equal or different}
  \item{\dots}{further arguments passed on to function \code{Mclust}}
}
\details{
  This function uses the package \code{mclust} to fit a mixture of two
  Gaussians to the data.
  The threshold is then estimated from the fitted Gaussian with the
  lower mean value.
}
\value{
  Single numeric value. The threshold that is the minimal
  data value with an adjusted p-value smaller than a specified
  threshold.
}
\author{Joern Toedling, Aleksandra Pekowska}
\note{
  Thanks to Richard Bourgon for pointing out the necessity of providing
  this method as an alternative way of estimating the threshold.
}
\seealso{\code{mclust}, \code{\link[stats]{p.adjust}}}
\examples{
  exDir <- system.file("exData",package="Ringo")
  load(file.path(exDir,"exampleProbeAnno.rda"))
  load(file.path(exDir,"exampleX.rda"))
  smoothX <- computeRunningMedians(exampleX, probeAnno=exProbeAnno,
     modColumn = "Cy5", allChr = "9", winHalfSize = 400)

  ## compare the two different ways of estimating the threshold
  y0a <- apply(exprs(smoothX), 2, upperBoundNull)
  y0b <- apply(exprs(smoothX), 2, twoGaussiansNull)

  hist(exprs(smoothX)[,1], n=10, main=NA,
       xlab="Smoothed expression level [log2]")
  abline(v=c(y0a, y0b), col=c("blue","orange"), lwd=2)
  legend(x="topright", col=c("blue","orange"), lwd=2, 
         legend=c(expression(paste(y[0]," Non-parametric")),
                  expression(paste(y[0]," Gaussian"))))
}
\keyword{manip}
