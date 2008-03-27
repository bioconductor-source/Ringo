\name{getFeats}
\alias{getFeats}
\alias{getFeatures}
\title{Utility function to extract all features from a cherList}
\description{
  This is a small utility function for extracting all related features
  from a \emph{cherList}, a list of \emph{ChIP-enriched regions}.
}
\usage{
getFeats(cl)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cl}{object of class \code{cherList}, a list of cher objects}
}
%\details{}
\value{
  a character vector containing the names of all features that were
  associated to any ChIP-enriched region in the list before, using the
  function \code{\link{relateChers}}
}
\author{Joern Toedling}
\seealso{\code{\link{relateChers}},\code{\link{cher-class}}}
%\examples{}
\keyword{manip}
