\name{qop-class}
\Rdversion{1.1}
\docType{class}
\alias{qop-class}
\alias{plot,qop-method}
\alias{plot,qop,ANY-method}
\alias{show,qop-method}
\title{Class "qop" Quantiles Over Positions}
\description{
  Class for storing result objects of the function
  \code{quantilesOverPositions}
}
\section{Objects from the Class}{
  Objects can be created by calls to the function
  \code{quantilesOverPositions}.
}
\section{Slots}{
	 \describe{
    \item{\code{data}:}{Object of class \code{"list"} }
    \item{\code{genes}:}{Object of class \code{"character"} }
    \item{\code{positions}:}{Object of class \code{"numeric"} }
    \item{\code{samplenames}:}{Object of class \code{"character"} }
    \item{\code{quantiles}:}{Object of class \code{"numeric"} }
    \item{\code{mapping}:}{Object of class \code{"character"} }
  }
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "qop")}: plot the data}
    \item{show}{\code{signature(object = "qop")}: give a short
      description of the object}
	 }
}
\author{Joern Toedling}
%\note{ ~~further notes~~ }
\seealso{\code{\link{quantilesOverPositions}}}
\examples{
  ## see the man page of the function
  ## 'quantilesOverPositions'
}
\keyword{classes}
\keyword{internal}
