## (c) 2007  Joern Toedling

###############################################################################
##  script to visualize enrichment/Chers-to-transcript
##  (or other genomic features) relations using Rgraphviz
###############################################################################

library("Ringo")
library("Rgraphviz")

doFunctionTest <- FALSE
analyzeTfChers3 <- FALSE

# small function to capitalize 1st letter of each element of character vector
capit <- function(x){ # capitalize first letter
  stopifnot(is.character(x))
  sapply(strsplit(x, ""), function(y)
         if (length(y)==0) return("")
         else paste(toupper(y[1]),paste(tolower(y[-1]),collapse=""), sep=""))
}#capit

## this function converts a cherList into a two-column matrix from a graph
## can easily be constructed using the Rgraphviz function 'ftM2graphNEL'
cherList2AssignTable <- function(pl, target.names=c("upSymbol","downSymbol"),
                                 tableNames=c("antibody","target"), verbose=TRUE)
{
  stopifnot(is.list(pl), inherits(pl[[1]],"cher"),
            "antibody" %in% slotNames(pl[[1]]), length(tableNames)==2)
  if (!all(target.names %in% names(pl[[1]]@extras)))
    stop(paste("Some of the target names", paste(target.names, collapse=", "),
               "are not specified for the chers in the list.\n"))
  resTable <- matrix(NA, nrow=length(pl)*10, ncol=2)
  nAssigns <- 1
  for (i in 1:length(pl)){
    if (verbose && (i %% 1000 ==0)) cat(i,"... ")
    this.cher <- pl[[i]]
    ## put upstream and downstream possible effects together
    these.targets <- unique(unlist(this.cher@extras[target.names], use.names=FALSE))
    these.targets <- these.targets[these.targets!=""]
    if (length(these.targets)<1) next
    for (this.target in these.targets){
      resTable[nAssigns,1:2] <- capit(c(this.cher@antibody, this.target))
      nAssigns <- nAssigns + 1
    }#for (this.sym)
  }# for i
  resTable <- resTable[complete.cases(resTable),,drop=FALSE]
  if (nrow(resTable)== 0)
    warning("No cher-target relations found.")
  resTable <- resTable[!duplicated(resTable),,drop=FALSE]
  resTable <- apply(resTable, 2 , function(x) gsub("[-_]",".",x))
  colnames(resTable) <- tableNames
  return(resTable)
}#cherList2AssignTable

## next function is a wrapper around tfM2graphNEL
assignTable2Graph <- function(assignTab, weights=NULL, antibody.colors=NULL, target.color="#B0C4DE", edge.color="darkslategrey", by.antibody=TRUE, sel.patterns=NULL, weight.colors=NULL, weight.breaks=NULL){
  # by.antibody = should edges by color-coded by used antibody?
  # sel.nodes = if specified, only a sub-graph is created in which each node names matches to any of these regular expression patterns (see ?grep)
  if (!is.null(weights))
    stopifnot(length(weights)==nrow(assignTab))
  myGraph <- ftM2graphNEL(assignTab, W=weights)
  if (!is.null(sel.patterns)){
    sel.nodes <- unique(unlist(lapply(as.list(sel.patterns), grep, x=nodes(myGraph), ignore.case=TRUE, value=TRUE)))
    myGraph <- subGraph(sel.nodes, myGraph)
  }# if (is.null(sel.patterns))
  ## remainder of function is setting attributes/colors for plotting of the graph
  my.antibodies <- unique(assignTab[,1])
  if (is.null(antibody.colors)){
    antibody.colors <- colorRampPalette(brewer.pal(9,"Spectral"))(length(my.antibodies))
    names(antibody.colors) <- my.antibodies
  } else if (length(antibody.colors)==1) {
    antibody.colors <- rep(antibody.colors, length(my.antibodies))
    names(antibody.colors) <- my.antibodies
  } else {
    stopifnot(all(my.antibodies %in% names(antibody.colors)))
    antibody.colors <- antibody.colors[my.antibodies]
  }
  ## set some graph attributes for coloring
  myGraphAttrs <- list(node = list(fillcolor=target.color, fontsize="20", fixedsize="false"), edge = list(color = edge.color))
  myNodeAttrs <- list()
  myNodeAttrs$fillcolor <- as.list(antibody.colors)
  myEdges <- edges(myGraph)
  myEdgeAttrs <- list()
  if (by.antibody|!is.null(weight.colors)){
    myEdgeAttrs$color <- list()
    if (!is.null(weight.colors)){
      if (is.null(weight.breaks))
        weight.breaks <- seq(min(weights),max(weights)+0.1,length.out=length(weight.colors)+1)
      else
        stopifnot(length(weight.breaks)==length(weight.colors)+1,
                  min(weights)>=min(weight.breaks), max(weights)<max(weight.breaks))
      myEdgeWeights <- edgeWeights(myGraph)
    } #if (!is.null(weight.colors))
    for (this.antibody in my.antibodies)
      for (this.target in myEdges[[this.antibody]])
        myEdgeAttrs$color[[paste(this.antibody,this.target,sep="~")]] <-
          ifelse(is.null(weight.colors),
                 antibody.colors[this.antibody],
                 weight.colors[findInterval(myEdgeWeights[[this.antibody]][this.target], weight.breaks)])
  }#if (by.antibody|!is.null(weight.colors))
  ### color edges by their weight:
  attr(myGraph,"graphAttrs") <- myGraphAttrs
  attr(myGraph,"nodeAttrs") <- myNodeAttrs
  attr(myGraph,"edgeAttrs") <- myEdgeAttrs
  #class(myGraph) <- c("cherAssignGraph",class(myGraph))
  return(myGraph)
}# assignTable2Graph

myGraphPlot <- function(x, edge.lwd=2, ...){
  stopifnot(inherits(x,"graphNEL"))
  par(mar=c(0.1,0.1,0.1,0.1), font=2, lwd=edge.lwd)
  plot(x, attrs=attr(x,"graphAttrs"), nodeAttrs=attr(x,"nodeAttrs"),
       edgeAttrs=attr(x,"edgeAttrs"),...)
}#myGraphPlot

### do testing of functions? might become example code in distant future
if (doFunctionTest) {
  example(findChersOnSmoothed)
  assignTabX <- cherList2AssignTable(chersX,target.name=c("upSymbol","downSymbol"))
  print(assignTabX)
  ### convert table into a graph without weights:
  graphX <- assignTable2Graph(assignTabX)
  myGraphPlot(graphX)
  ### convert table into a graph with weights:
  graphXw <- assignTable2Graph(assignTabX, weights=c(1,1,-1,-1), weight.colors=c("green","red"), weight.breaks=c(-1,0,1.1))
  myGraphPlot(graphXw)
}#if (doFunctionTest)


#############################################################################
### the following is how we used these functions with our TF ChIP-chip data,
### feel free to modify it to suit your own data
#############################################################################

if (analyzeTfChers3){
### tell R in which directory the saved binary data resides
dataDir <- "/panfs/panda1/huber/sperling/tfchip/data"

### load list of identified enrichments:
load(file=file.path(dataDir,"chersV3c.RData"))
## the loaded object allTFChers is a list of lists of chers
tfChers <- unlist(allTFChers,recursive=FALSE, use.names=FALSE)
if ("modification" %in% names(tfChers[[1]]))
  tfChers <- lapply(tfChers, function(p){names(p)[names(p)=="modification"] <- "antibody"; p})

assignTF <- cherList2AssignTable(tfChers)
assignTF[,1] <- gsub("H1.cdna.2","Dpf3.1", assignTF[,1])
assignTF[,1] <- gsub("H1.cdna.3","Dpf3.2", assignTF[,1])

my.patterns <- unique(c("Nkx2","Myocd", "GATA4","SRF$", "MEF2","bcl2$","Hand","Alk2","Irx4","Tbx","Isl1","Bop","Anf","Dpf3$", paste("^",names(which(table(assignTF[,2])>5)),"$",sep="")))
tf.colors <- brewer.pal(8,"Dark2"); tf.colors[4] <- "#FF00FF"; tf.colors[3] <- "royalblue";
names(tf.colors) <- c("Gata4","Nkx2.5","Mef2","Srf","","Dpf3.1","Dpf3.2","Anti.flag")
graphTF <- assignTable2Graph(assignTF, sel.patterns=my.patterns, antibody.colors=tf.colors)

#pdf(file.path(dataDir,"graphTF.pdf"),width=12,height=10)
myGraphPlot(graphTF)
#dev.off()
}#if (analyzeTfChers3)
