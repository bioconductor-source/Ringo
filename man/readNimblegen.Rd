\name{readNimblegen}
\alias{readNimblegen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function to read in Nimblegen Intensity Text Files}
\description{
  Function to read in Nimblegen Intensity Text Files into an
  \code{RGList}. Calls some other functions for actual reading of data.
  This function is to be used with two-color NimbleGen array data.
  Use the function \code{read.xysfiles} of the \code{oligo} package
  for single-color data.
}
\usage{
readNimblegen(hybesFile, spotTypesFile, path = getwd(),
 headerPattern="# software=NimbleScan",verbose = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hybesFile}{Name of the file describing the arrays. In
    \code{limma} this file would be called targets file.}
  \item{spotTypesFile}{spot types also used by \code{limma}}
  \item{path}{Path to directoy that hold the files \code{hybesFile},
    \code{spotTypesFile} and also the intensity files. Set this to
    \code{NULL} if you prefer the arguments \code{hybesFile},
    \code{spotTypesFile} and the file-name entries of the hybes file to 
    be treated as absolute or relative file paths themselves.}
  \item{headerPattern}{string; pattern used to identify explantory
    header lines in the supplied pair-format files}
  \item{verbose}{logical; progress output to STDOUT?}
  \item{\dots}{further arguments passed on the \code{readNgIntensitiesTxt}}
}
\value{
  Returns raw intensity values in form of an \code{RGList}.
}
\author{Joern Toedling \email{toedling@ebi.ac.uk}}
\seealso{\code{\link[limma]{rglist}},
  \code{\link[limma]{readTargets}}}
\examples{
  exDir <- system.file("exData",package="Ringo")
  exRG <- readNimblegen("example_targets.txt","spottypes.txt",path=exDir)
  print(exRG)
}
\keyword{file}% at least one, from doc/KEYWORDS
\keyword{IO}% at least one, from doc/KEYWORDS
