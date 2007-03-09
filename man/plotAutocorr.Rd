\name{plot.autocor.result}
\alias{plotAutocor}
\alias{plot.autocor.result}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plots auto-correlation of probe intensities}
\description{
  Function to plot the auto-correlation of probe intensities
  computed by function \code{autocor}.
}
\usage{
\method{plot}{autocor.result}(x, 
plot.title = "ChIP: Autocorrelation of Intensities", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of class \code{autocor.result}, the result  of
    function \code{autocor}}
  \item{plot.title}{\code{main} title of the plot}
  \item{\dots}{further arguments passed on to \code{plot.default}, see details}
}
\details{
  The following arguments to \code{plot.default} are already defined in
  the function and thus cannot be specified by the user as further
  arguments in \dots:
  \code{type, lwd, xlab, ylab, col}.
  Argument \code{main} is specified in \code{plot.title}. 
}
\value{
  invisible \code{NULL}
}
\author{Joern Toedling \email{toedling@ebi.ac.uk}}
\seealso{\code{\link{autocor}}}
\examples{
  ## see the help page of 'autocor' for an example
}
\keyword{hplot}
