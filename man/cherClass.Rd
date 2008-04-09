\name{cher-class}
\docType{class}
\alias{cher-class}
\alias{initialize,cher-method}
\alias{show,cher-method}
\alias{update,cher-method}
\alias{cher}
\alias{Cher}
\alias{cherList}
\alias{cherList-class}
\title{Class "cher" - ChIP-enriched region}
\description{
  An object of class \code{cher} (ChIP-enriched region) holds
  characteristics of an enriched genomic region from ChIP chip
  data. }
\section{Objects from the Class}{
  Objects can be created by calls of the form
  \code{new("cher", name, chromosome, start, end, cellType, antibody,
    maxLevel, score, probes, extras, ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{name}:}{character vector of length 1 unequivocally
      describing the cher, e.g. "Suz12.Nudt2.upstream.cher"}
    \item{\code{chromosome}:}{character vector of length one, naming the
      chromosome of the region, e.g. "9"}
    \item{\code{start}:}{\code{integer}, region start position on the
      chromosome, e.g. 34318900} 
    \item{\code{end}:}{\code{integer}, region end position on the
      chromosome, e.g. 34320100}
    \item{\code{cellType}:}{\code{character} vector describing the cell
      type the ChIP chip experiment has been done on,
      e.g. "HeLa" or "human"}
    \item{\code{antibody}:}{\code{character} vector describing the
      antibody or characteristic for which fragments were supposedly
      enriched in immuno-precipitation step,
      e.g. "Suz12" for the protein Suz12}
    \item{\code{maxLevel}:}{\code{numeric}, maximal (smoothed) probe
      level in the cher, e.g. 2.00}
    \item{\code{score}:}{\code{numeric} of a cher score, currently we
      use the sum of smoothed probe levels (log fold changes),
      e.g. 69.16}
    \item{\code{probes}:}{\code{character vector} of probe identifiers
      of all probes with match positions in the \code{cher}}
    \item{\code{extras}:}{\code{list} of further elements used to
      annotate the \code{cher}; examples of such that are used in
      \code{Ringo} are:
      \describe{
	\item{typeUpstream}{optional character vector of features that
	  this \code{cher} is located upstream of,
	  e.g. the transcriptional start site of "ENST00000379158".
	  See \code{\link{relateChers}} for details.}
	\item{typeInside}{ optional character vector of features
	  that this \code{cher} is located inside of }
	\item{distMid2TSS}{optional named numeric vector of
	  distances of the \code{cher}'s middle position to
	  features, e.g. TSSs of features upstream and inside; names
	  are the features to which the distances are given; only
	  meaningful in combination with \code{typeUpstream} and
	  \code{typeInside}; e.g. 55 with name "ENST00000379158"}
	\item{upSymbol}{optional character vector of gene symbols of
	  features the cher is located upstream of; supplements
	  \code{typeUpstream}; e.g. "Nudt2" }
	\item{inSymbol}{optional character vector of gene symbols of
	  features the cher is located upstream of; supplements
	  \code{typeInside}.}
	\item{\ldots}{further list elements can be added using the
	  \code{update} method.}
      }
    }
  }
}
\section{Methods}{
  \describe{
    \item{initialize}{create a new \code{cher}; see section
      \code{examples} below}
    \item{plot}{calls \code{\link{chipAlongChrom}} to plot the
      \code{cher}; see \code{\link{plot.cher}} for more details}
    \item{update}{signature(\code{cher},...); updates elements of the
      \code{cher} object; The further arguments in '...' are
      interpreted. Arguments corresponding to defined slot names of the
      \code{cher} result in the value by that slot being replaced by the
      specified value for the argument; Argment names that do not
      correspond to slot names of the object result in list elements of
      the \code{extras} list of the \code{cher} being replaced by the
      given values for these arguments or the values are appended to the
      current \code{extras} list and the argument names make up the list
      names of the appended arguments. See section
      \code{examples} below for an example how to use this method.}
  }
}
\author{Joern Toedling, Tammo Krueger}
\note{
  The \code{cher} class used to be an S3 list before.

  The term 'cher' is shorthand for 'ChIP-enriched region'.
  We think this term is more appropriate than the term
  'peak' commonly used in ChIP-chip context. Within such regions the
  actual signal could show two or more actual signal peaks or none at
  all (long plateau). 
}
\section{cherList}{
  A list in which each element is of class \code{cher}, is called
  a \code{cherList}. This class, however, is rarely used (yet).
}
\seealso{\code{\link{plot.cher}}, \code{\link{findChersOnSmoothed}}, \code{\link{relateChers}}}
\examples{
  ## how to create a cher object from scratch
  cherNudt2 <- new("cher", name="nudt2.cher", chromosome=9,
                   start=34318954, end=34319944, antibody="Suz12",
                   maxLevel=2.00, score=69.2, upSymbol="NUDT2")
                   #extras=list(upSymbol="NUDT2"))
  cherNudt2
  str(cherNudt2)

  ## use the update method (note:this update is biologically meaningless)
  cher2 <- update(cherNudt2, cellType="HeLa", downSymbol="P53",
                  probes=c("probe1","probe2"))
  cher2; str(cher2)

  ## plot a cher object
  exDir <- system.file("exData",package="Ringo")
  load(file.path(exDir,"exampleProbeAnno.rda"))
  load(file.path(exDir,"exampleX.rda"))
  smoothX <- computeRunningMedians(exampleX, probeAnno=exProbeAnno,
       modColumn = "Cy5", allChr = "9", winHalfSize = 400)
  plot(cherNudt2, smoothX, probeAnno=exProbeAnno, gff=exGFF, extent=5000)
}
\keyword{classes}
