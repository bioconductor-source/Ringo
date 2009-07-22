### function for using topGO with a supplied list of genes:
sigGOTable <- function(selGenes, gene2GO, universeGenes, ontology="BP", maxP=0.001)
{
  require("topGO") # for new version: >= 1.13.1
  stopifnot(is.character(selGenes), is.numeric(maxP))
  ontology <- match.arg(ontology, c("BP","MF","CC"))
  if (missing(universeGenes))
    universeGenes <- names(gene2GO)
  inGenes <- factor(as.integer(universeGenes %in% selGenes))
  names(inGenes) <- universeGenes
  GOdata <- new("topGOdata", ontology=ontology, allGenes=inGenes, 
                annot=annFUN.gene2GO, gene2GO=gene2GO)
  resultFisher <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
  sTab <- GenTable(GOdata, p.value=resultFisher, topNodes=length(usedGO(GOdata)))
  sTab <- subset(sTab, as.numeric(p.value) < maxP)
  return(sTab)
}
