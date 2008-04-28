\name{validProbeAnno}
\alias{validProbeAnno}
\alias{validateProbeAnno}
\title{Function to check a probeAnno environment}
\description{
  This function checks whether a probeAnno environment seems to be valid
  and will work with other Ringo functions.
}
\usage{
validProbeAnno(probeAnno)
}
\arguments{
  \item{probeAnno}{an environment that holds probe matches to genomic
    coordinates}
}
\details{
  This function checks certain properties of the mapping environment
  that are used by other Ringo functions. It can indicate potential
  problems in the environment.
}
\value{
  \code{TRUE} if the environment seems to be a valid probeAnno
  environment, \code{FALSE} if a potential problems with this
  environment was identified. This potential problem is explained as a
  warning.
}
\author{Joern Toedling \email{toedling@ebi.ac.uk}}
\seealso{\code{\link{posToProbeAnno}}}
\examples{
 ## first a toy example:
 if (interactive()){
  A = new.env()
  assign("1.+.start",seq(100,1000,by=100),env=A)
  validProbeAnno(A)
  assign("1.+.end",c(99,seq(250,1050,by=100)),env=A)
  assign("1.+.unique",numeric(10),env=A)
  assign("1.+.index",c(2:5,1,7,8,6,10), env=A)
  validProbeAnno(A)
  assign("1.+.index",c(2:5,1,7,8,6,10,3), env=A)
  validProbeAnno(A)
  assign("1.+.end",c(150,seq(250,1050,by=100)),env=A)
  validProbeAnno(A)
 } 
 ## validate the provided example probeAnno
 load(file.path(system.file("exData",package="Ringo"),"exampleProbeAnno.rda"))
 validProbeAnno(exProbeAnno)
}
\keyword{internal}
