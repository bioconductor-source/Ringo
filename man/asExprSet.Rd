\name{asExprSet}
\alias{asExprSet}
\alias{asExpressionSet}
\title{converts a Ringo MAList into an ExpressionSet}
\description{
  Function to convert an object of class \code{MAList} into an object of
  class \code{ExpressionSet}. Note that the otherwise optional
  \code{targets} component is required in this case to generate the
  \code{phenoData} of the new \code{ExpressionSet}.
}
\usage{asExprSet(from, idColumn="PROBE_ID")}
\arguments{
  \item{from}{object of class \code{MAList} to convert into an
    \code{ExpressionSet}}
  \item{idColumn}{string; indicating which column of the \code{genes}
    data.frame of the MAList holds the identifier for reporters on the
    microarray. This column, after calling
    \code{\link[base]{make.names}} on it, will make up the unique
    \code{featureNames} of the resulting \code{ExpressionSet}.}
}
\value{
  an object of class \code{ExpressionSet}
}
\note{
  There is a more general function for converting MALists to
  ExpressionSets in the package \code{convert}. This function here is
  solely intended for converting Ringo-generated \code{MAList}s into
  \code{ExpressionSet}s.
}
\author{Joern Toedling}
\seealso{\code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}, \code{\link{preprocess}}}
\examples{
  exDir <- system.file("exData",package="Ringo")
  exRG  <- readNimblegen("example_targets.txt","spottypes.txt",path=exDir)
  exMA  <- preprocess(exRG, "none", returnMAList=TRUE)
  exX   <- asExprSet(exMA)
}
\keyword{manip}% at least one, from doc/KEYWORDS
