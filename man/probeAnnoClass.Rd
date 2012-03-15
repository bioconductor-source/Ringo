\name{probeAnno-class}
\docType{class}
\alias{probeAnno-class}
\alias{probeAnno}
\alias{arrayName<-,probeAnno,character-method}
\alias{arrayName,probeAnno-method}
\alias{arrayName<-}
\alias{arrayName}
\alias{genome<-}
\alias{genome}
\alias{chromosomeNames}
\alias{[<-,probeAnno-method}
\alias{[,probeAnno-method}
\alias{[,probeAnno,ANY,ANY,ANY-method}
\alias{[<-,probeAnno,ANY,ANY,ANY-method}
\alias{chromosomeNames,probeAnno-method}
\alias{get,character,missing,probeAnno-method}
\alias{initialize,probeAnno-method}
\alias{ls,probeAnno,missing,missing,missing,missing-method}
\alias{show,probeAnno-method}
\alias{genome<-,probeAnno,character-method}
\alias{genome,probeAnno-method}
\alias{coerce,environment,probeAnno-method}
\title{Class "probeAnno"}
\description{A class that holds the mapping between reporters/probes
  on a microarray and their genomic match  position(s) in a chosen
  genome.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("probeAnno",
    map, arrayName, genome)}.
}
\section{Slots}{
  \describe{
    \item{\code{map}:}{Object of class \code{"environment"} This map
      consists of four vectors for each chromosome/strand, namely, say
      for chromosome 1:
      \describe{
	\item{1.start}{genomic start coordinates of probe matches on
	  chromosome 1}
	\item{1.end}{genomic start coordinates of probe matches on
	  chromosome 1}
	\item{1.index}{identifier of probes matching at these
	  coordinates}
	\item{1.unique}{vector of the same length as the three before;
	  encoding how many matches the corresponding probe has in the
	  given file or \code{data.frame}. An entry of '0' indicates
	  that the probe matching at this position has only this one
	  match.}
    }}
    \item{\code{arrayName}:}{Object of class \code{"character"}, the
      name or identifier of the microarray design,
      e.g. \code{2005-06-17\_Ren\_MM5Tiling\_Set1}}
    \item{\code{genome}:}{Object of class \code{"character"}, which
      genome the reporters have been mapped to}
  }
}
\section{Methods}{
  \describe{
    \item{arrayName}{obtain the microarray platform name}
    \item{arrayName<-}{replace the microarray platform name}
    \item{[}{get elements from the \code{map} environment}
    \item{[<-}{assign elements to the \code{map} environment}
    \item{chromosomeNames}{obtain a character vector holding the names
      of the chromosomes for which the probeAnno objects holds a mapping.}
    \item{get}{get elements from the \code{map} environment}
    \item{initialize}{create mew probeAnno object}
    \item{ls}{list elements of the \code{map} environment }
    \item{genome}{obtain the description of the genome the reporters
      were mapped to}
    \item{genome<-}{replace the description of the genome the reporters
      were mapped to}
    \item{as}{signature(from="environment"); function to coerce
      old-style 'probeAnno' environments to new-style 'probeAnno'
      objects. Simply creates a new object with the old environment in
      its \code{map} slot}
  }
}
\author{Joern Toedling; Wolfgang Huber}
\note{'probeAnno' objects used to be environments and still are used as
  such in package \code{tilingArray}}
\seealso{posToProbeAnno}
\examples{
  pa <- new("probeAnno")
  pa["X.start"] <- seq(5000,10000,by=1000)
  if (interactive()) show(pa)
  pa2 <- posToProbeAnno(file.path(system.file("exData",package="Ringo"),
                                  "MOD_2003-12-05_SUZ12_1in2.pos"))
  arrayName(pa2) <- "NimbleGen MOD_2003-12-05_SUZ12_1in2"
  genome(pa2) <- "H. sapiens (hg18)"
  show(pa2)
  head(pa2["9.start"])
}
\keyword{classes}
