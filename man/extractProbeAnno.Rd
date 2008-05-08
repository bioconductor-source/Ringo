\name{extractProbeAnno}
\alias{extractProbeAnno}
\alias{probeAnnoFromRGList}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Build probeAnno from match positions in an RGList}
\description{
  This function can be used to build a probeAnno object from the
  reporter match positions given in the  'genes' slot of an RGList if
  present, as is the case with some ChIP-chip microarray platforms,
  e.g. with certain Agilent ones after reading in the data with
  \code{read.maimages(...,"agilent")}.
}
\usage{
extractProbeAnno(object, format = "agilent", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object that holds the data and the probe match
    positions, currently can only be of class \code{RGList}}
  \item{format}{in which format are the reporter match positions stored
    in the object; see details; currently only "agilent" is implemented}
  \item{\dots}{further arguments that are passed on to the function
    \code{posToProbeAnno}}
}
\details{
  Which information is used for creating the probeAnno is specified by
  the argument \code{format}.
  \describe{
    \item{agilent}{
      expects that the object is of class \code{RGList}. The 'genes'
      element of the object is taken. This element is expected to have
      at least a column 'ProbeName', which stores the unique
      reporter/probe identifiers, and a column 'SystematicName', which
      holds the probe match position in the format
      "chr<chromosome>:coordinate1-coordinate2",
      e.g."chr1:087354051-087354110".
    }
  }
}
\value{
  An object of  class \code{probeAnno} holding the mapping between
  reporters and genomic positions.
}
\author{Joern Toedling \email{toedling@ebi.ac.uk}}
\seealso{ \code{\link{posToProbeAnno}}, \code{\link{probeAnno-class}}}
%\examples{}
\keyword{manip}
