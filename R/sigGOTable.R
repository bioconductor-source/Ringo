### function for using topGO with a supplied list of genes:
sigGOTable <- function(selGenes, universeGenes, annotType="list",
                       ontology="BP", maxP=0.001, algorithm="elim",
                       minGenes=10, gene2GO=list(), pkg="org.Mm.eg.db", ...)
{
  require("topGO") # version: >= 1.15.0
  stopifnot(is.character(selGenes), is.numeric(maxP))
  ontology  <- match.arg(ontology, c("BP","MF","CC"))
  annotType <- match.arg(annotType, c("list", "org"))
  if (annotType=="org")
    require(pkg, character.only=TRUE)
  algorithm <- match.arg(algorithm,
                         c("elim","classic","weight","topgo","parentChild"))
  if (missing(universeGenes)){
    if (annotType=="list")
      universeGenes <- names(gene2GO)
    else if (annotType=="org")
      universeGenes <- mappedkeys(
         get(paste(gsub(".db$","", pkg), "GO", sep="")) )
  }
  inGenes <- factor(as.integer(universeGenes %in% selGenes))
  names(inGenes) <- universeGenes
  if (annotType=="list")
    GOdata <- new("topGOdata", ontology=ontology,
                  allGenes=inGenes, annot=annFUN.gene2GO,
                  gene2GO=gene2GO, nodeSize=minGenes, ...)
  if (annotType=="org")                 
    GOdata <- new("topGOdata", ontology=ontology,
                  allGenes=inGenes, nodeSize=minGenes,
                  annot=annFUN.org, mapping=pkg, ...)
  resultFisher <- runTest(GOdata, algorithm=algorithm,
                          statistic="fisher", cutOff=maxP)
  sTab <- GenTable(GOdata, p.value=resultFisher,
                   topNodes=length(usedGO(GOdata)))
  sTab <- subset(sTab, as.numeric(p.value) < maxP)
  return(sTab)
}
