\name{asExprSet}
\alias{asExprSet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{converts MAList into an exprSet}
\description{
  Function to convert an object of class \code{MAList} into an object of
  class \code{ExpressionSet}. Note that the optional \code{targets} component
  is needed in this case to generate the \code{phenoData} of the new
  \code{ExpressionSet}.
}
\usage{asExprSet(myMA)}
\arguments{
  \item{myMA}{object of class \code{MAList}}
}
\value{
  an object of class \code{ExpressionSet}
}
\author{Joern Toedling}
\seealso{\code{\link[Biobase]{ExpressionSet}}}
\examples{
  exDir <- system.file("exData",package="Ringo")
  exRG  <- readNimblegen("example_files.txt","spottypes.txt",path=exDir)
  exMA  <- preprocess(exRG, "none", returnMAList=TRUE)
  exX   <- Ringo:::asExprSet(exMA)
}
\keyword{internal}% at least one, from doc/KEYWORDS
