\name{ftr2xys}
\alias{ftr2xys}
\title{Convert a NimbleScan ftr-file into a xys-file}
\description{
  Auxiliary function to convert a NimbleScan \code{feature-report
    file} into a xys-file that can be used with the function
  \code{read.xysfiles} of package \code{oligo}.
}
\usage{
ftr2xys(ftr.file, path=getwd())
}
\arguments{
  \item{ftr.file}{character; file path of feature report file to convert
    into an xys file}
  \item{path}{file path to directory where the xys-file should be
    written to; defaults to the current working directory}
}
\details{
  The output file is names as the input ftr file; with the file
  extension \code{.ftr} replaced by \code{.xys}.
}
\value{
  Function returns only \code{NULL} invisibly and is only called for
  its side effect to write the xys-file into the current
  working directory.
}
\author{Joern Toedling}
\note{
  This function should only be used with one-color Nimblegen
  microarrays and when the correct xys-file of the raw
  data is not available. The output file can be used with the function
  \code{read.xysfiles} of package \code{oligo}.
}
\examples{
 \dontrun{
   sapply(list.files(pattern=".ftr$"),ftr2xys)
   library(oligo)
   fs = read.xysfiles(list.xysfiles())
 }
}
\keyword{file}
\keyword{IO}
