\name{exportCherList}
\alias{exportCherList}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function to export cherList into a file}
\description{
 Function to export cherList into a file of gff or BED format. This
 files can be imported as tracks into genome browsers.
}
\usage{
exportCherList(object, filename = "chers.gff", format = "gff3", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object of class \code{cherList}}
  \item{filename}{character; path to file to be written}
  \item{format}{Format of exported file; currently only "gff3" and "bed"
    are supported}
  \item{\dots}{further arguments to be passed on to the \code{trackSet}
    method}
}
\details{
  First converts the \code{cherList} into an object of class
  \code{trackSet} from package \pkg{rtracklayer} and then calls the
  \code{export} method as defined for a \code{trackSet}.
}
\value{
  returns invisible NULL; called for the side effect of writing the file
  \code{filename}.
}
\author{Joern Toedling}
\seealso{Class \code{trackset} in package \pkg{rtracklayer}}
\examples{
\dontrun{
  exDir <- system.file("exData",package="Ringo")
  load(file.path(exDir,"exampleProbeAnno.rda"))
  load(file.path(exDir,"exampleX.rda"))
  smoothX <- computeRunningMedians(exampleX, probeAnno=exProbeAnno,
       modColumn = "Cy5", allChr = "9", winHalfSize = 400)
  chersX <- findChersOnSmoothed(smoothX, probeAnno=exProbeAnno,
       thresholds=0.45, allChr="9", distCutOff=600, cellType="human")
  exportCherList(chersX, file="chers.gff")
}
}
\keyword{IO}
\keyword{file}
% __ONLY ONE__ keyword per line
