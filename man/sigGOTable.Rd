\name{sigGOTable}
\Rdversion{1.1}
\alias{sigGOTable}
\title{
  Obtain significant GO terms for a list of genes
}
\description{
  Function to obtain a table of significant GO terms given a list of
  genes, for example: genes related to ChIP-enriched regions.
}
\usage{
sigGOTable(selGenes, universeGenes, annotType="list",
           ontology="BP", maxP=0.001, algorithm="elim",
           minGenes=10, gene2GO=list(), pkg="org.Mm.eg.db", ...)
}
\arguments{
  \item{selGenes}{character; names of the genes in the list; must
    correspond to the names of the list \code{gene2GO} (if annotType is
    \dQuote{list}) or to Entrezgene identifiers (if annotType is
    \dQuote{org}). }
  \item{universeGenes}{
    character; names of the genes to focus the analysis on (if
    desired). Sometimes specifying a smaller list of genes as the
    background 'universe' may lead to more informative results from the
    Fisher test. }
  \item{annotType}{string; in which format is the GO annotation of the
    genes provided; either \dQuote{list} if a named list (argument
    \code{gene2GO}) is provided OR \dQuote{org} if one of the
    Entrezgene-based annotation packages, e.g. \code{org.Xx.eg.db},
    shall be used }
  \item{ontology}{string; which ontology to use, one of \dQuote{BP},
    \dQuote{MF}, or \dQuote{CC} }
  \item{maxP}{numeric; which maximum p-value to allow for GO terms to be
    called significantly associated to the gene list; passed on topGO
    function \code{runTest} and used to subset the table of significant
    GO terms reported by topGO. }
  \item{algorithm}{string; one of \sQuote{elim},  \sQuote{classic},
    \sQuote{weight}, \sQuote{topgo}, or \sQuote{parentChild}; which
    method to use in the topGO analyis; see the topGO package for
    details. }
 \item{minGenes}{numeric; the required minimum number of genes
   (from \code{universeGenes}) that must be annotated for a GO term in
   order to consider it in the analysis; corresponds to the argument
   \code{nodeSize} in \code{topGO} }
  \item{gene2GO}{
    A list that contains for each gene the identifiers (GO:00001234...)
    of GO terms that are annotated for this gene; only used if
    \code{annotType=}\dQuote{list} }
 \item{pkg}{Name of the Entrezgene-based annotation package to obtain
   the GO annotation from. Only used if \code{annotType=}\dQuote{org} }
 \item{\dots}{additional arguments that are passed on to the constructor
   when building the \code{topGOdata} object.}
}
\details{
  This function is really just a convenience wrapper around functions in
  the \code{topGO} package.

  It performs an exact Fisher-test and by default uses the \sQuote{elim}
  method of the \code{topGO} testing framework.
}
\value{
  a \code{data.frame} containing the GO terms with p-values lower than
  the specified cut-off.
}
\author{Joern Toedling}
\seealso{function \code{runTest} and the vignette in package
  \code{topGO}}
\examples{
 \dontrun{
   library("topGO")
   library("hgu95av2.db")
   # get probe sets annotated for 'regulation of blood pressure'
   sel <- get("GO:0007616", hgu95av2GO2ALLPROBES)
   # create list that matches probe sets to GO
   ps2GO <- lapply(as.list(hgu95av2GO), names)

   sigGOTable(selGenes=sel, annotType="list", gene2GO=ps2GO)
 }
}
\keyword{internal}
\keyword{manip}
