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
sigGOTable(selGenes, gene2GO, universeGenes, ontology = "BP", maxP = 0.001)
}
\arguments{
  \item{selGenes}{character; names of the genes in the list; must
    correspond to names of the list \code{gene2GO}
  }
  \item{gene2GO}{
    A list that contains for each gene the identifiers (GO:00001234...)
    of GO terms that are annotated for this gene.
  }
  \item{universeGenes}{
    character; names of the genes to focus the analysis on (if
    desired). Sometimes specifying a smaller list of genes as the
    background 'universe' may lead to more informative results from the
    Fisher test.
  }
  \item{ontology}{string; which ontology to use, one of \dQuote{BP},
    \dQuote{MF}, or \dQuote{CC}.
  }
  \item{maxP}{numeric; which maximum p-value to allow for GO terms to be
    called significantly associated to the gene list
  }
}
\details{
  This function is really just a convenience wrapper around functions in
  the \code{topGO} package.

  It performs an exact Fisher-test and uses the \sQuote{elim} method of
  the \code{topGO} testing framework.
}
\value{
  a \code{data.frame} containing the GO terms with p-values lower than
  the specified cut-off.
}
\author{Joern Toedling}
\note{
  This function in this form only works with the package
  topGO versions >= 1.13.1.
}
%\seealso{}
\examples{
 \dontrun{
   library("topGO")
   library("hgu95av2.db")
   # get probe sets annotated for 'regulation of blood pressure'
   sel <- get("GO:0007616", hgu95av2GO2ALLPROBES)
   # create list that matches probe sets to GO
   ps2GO <- lapply(as.list(hgu95av2GO), names)

   sigGOTable(selGenes=sel, gene2GO=ps2GO)
 }
}
\keyword{internal}
\keyword{manip}
