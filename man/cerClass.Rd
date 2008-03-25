\name{cher}
\alias{Cher}
\alias{cher-class}
\alias{print.cher}
\alias{cherList-class}
\alias{cherList}
\alias{ChIP-enriched region}
\docType{class}
\title{Class cher - A list holding ChIP-enriched region characteristics}
\description{
  An object of class \code{cher} is a list holding characteristics of an
  enriched genomic region from ChIP chip data. While this is only a S3
  rather than a formal \code{S4} class, the presence of certain list
  elements enables a standardized approach to working with enriched
  regions identified from ChIP-chip data.
  The term "cher" is shorthand for "ChIP-enriched region".
  We think this term is more appropriate than the term
  'peak' commonly used in ChIP-chip context. Within such regions the
  actual signal could show two or more actual signal peaks or none at
  all (long plateau). 
}
\section{Slots/List Components}{
  \describe{
    \item{name}{required character vector of length 1
    unequivocally describing the cher, e.g. "Suz12.Nudt2.upstream.cher"}
    \item{chr}{required character vector of length one, naming the
    chromosome of the region, e.g. "9" }
    \item{start}{required integer of region start position on the
    chromosome, e.g. 34318900 }
    \item{end}{required integer of region start position on the
    chromosome, e.g. 34320100 }
    \item{cellType}{optional character vector describing the cell
    type the ChIP chip experiment has been done on, e.g. "human"}
    \item{antibody}{required character vector describing the
    antibody or characteristic for which fragments were supposedly
    enriched in immuno-precipitation step,
    e.g. "Suz12" for the protein Suz12}
    \item{typeUpstream}{optional character vector of features that this
    cher is located upstream of,
    e.g. the transcriptional start site of "ENST00000379158".
    See \code{\link{relateChers}} for details.}
    \item{typeInside}{ optional character vector of features
    that this cher is located inside of }
    \item{maxCher}{optional numeric of maximal (smoothed) probe level
    in the cher, e.g. 2.00}
    \item{score}{optional numeric of a cher score, currently we use the
    sum of smoothed probe levels (log fold changes), e.g. 69.16}
    \item{probes}{optional character vector of probe identifiers of
    all probes with match positions in the cher }
    \item{distMid2TSS}{optional named numeric vector of distances of
    the cher middle position to features,
    e.g. TSSs of features upstream and inside; names are
    the features to which the distances are given; only
    meaningful in combination with \code{typeUpstream} and
    \code{typeInside}; e.g. 55 with name "ENST00000379158"}
    \item{upSymbol}{optional character vector of gene symbols of
    features the cher is located upstream of; supplements
    \code{typeUpstream}; e.g. "Nudt2" }
    \item{inSymbol}{optional character vector of gene symbols of
    features the cher is located upstream of; supplements
    \code{typeInside}.}
    \item{\ldots}{further list elements can be added like normal
      list elements.}
 }
}
\section{Creation}{
  Objects of class \code{cher} are created by other functions,
  such as \code{\link{findChersOnSmoothed}}; or by creating a list with
  the respective components and assigning it the class \code{cher},
  which is essentially, what the function \code{\link{newCher}} does.
}
\section{Methods}{
  Methods defined for lists of class \code{cher} are
  \code{print.cher} and \code{\link{plot.cher}}.
}
\section{cherList}{
  A list in which each element is of class \code{cher}, is called
  a \code{cherList}. This class, however, is not much in use yet.
}
\seealso{
  \code{\link{plot.cher}}, \code{\link{findChersOnSmoothed}},
  \code{\link{relateChers}}
}
\author{Joern Toedling, Tammo Krueger, PT Bote}
\examples{
  exDir <- system.file("exData",package="Ringo")
  load(file.path(exDir,"exampleProbeAnno.rda"))
  load(file.path(exDir,"exampleX.rda"))
  smoothX <- computeRunningMedians(exampleX, probeAnno=exProbeAnno,
   modColumn = "Cy5", allChr = "9", winHalfSize = 400)
  cersX <- findChersOnSmoothed(smoothX, probeAnno=exProbeAnno,
   thresholds=0.45, allChr="9", distCutOff=600)
  cersX <- relateChers(cersX, exGFF)
  str(cersX[[1]])
}
\keyword{classes}
\keyword{data}
