## (c) 2007  Joern Toedling

###############################################################################
##  script to visualize enrichment/peaks-to-transcript
##  (or other genomic features) relations using Rgraphviz
###############################################################################

library("Ringo")
library("Rgraphviz")

doFunctionTest <- FALSE

# small function to capitalize 1st letter of each element of character vector
capit <- function(x){ # capitalize first letter
  stopifnot(is.character(x))
  sapply(strsplit(x, ""), function(y)
         if (length(y)==0) return("")
         else paste(toupper(y[1]),paste(tolower(y[-1]),collapse=""), sep=""))
}#capit

## this function converts a peakList into a two-column matrix from a graph
## can easily be constructed using the Rgraphviz function 'ftM2graphNEL'
peakList2AssignTable <- function(pl, target.names=c("upSymbol","downSymbol"),
                                 tableNames=c("antibody","target"), verbose=TRUE)
{
  stopifnot(is.list(pl), inherits(pl[[1]],"peak"),
            "modification" %in% names(pl[[1]]), length(tableNames)==2)
  if (!all(target.names %in% names(pl[[1]])))
    stop(paste("Some of the target names", paste(target.names, collapse=", "),
               "are not specified for the peaks in the list.\n"))
  resTable <- matrix(NA, nrow=length(pl)*10, ncol=2)
  nAssigns <- 1
  for (i in 1:length(pl)){
    if (verbose && (i %% 1000 ==0)) cat(i,"... ")
    this.peak <- pl[[i]]
    ## put upstream and downstream possible effects together
    these.targets <- unique(unlist(this.peak[target.names], use.names=FALSE))
    these.targets <- these.targets[these.targets!=""]
    if (length(these.targets)<1) next
    for (this.target in these.targets){
      resTable[nAssigns,1:2] <- capit(c(this.peak$antibody, this.target))
      nAssigns <- nAssigns + 1
    }#for (this.sym)
  }# for i
  resTable <- resTable[complete.cases(resTable),,drop=FALSE]
  if (nrow(resTable)== 0)
    warning("No peak-target relations found.")
  resTable <- resTable[!duplicated(resTable),,drop=FALSE]
  resTable <- apply(resTable, 2 , function(x) gsub("[-_]",".",x))
  colnames(resTable) <- tableNames
  return(resTable)
}#peakList2AssignTable

## next function is a wrapper around tfM2graphNEL
assignTable2Graph <- function(assignTab, antibody.colors=NULL, target.color="#B0C4DE", edge.color="darkslategrey", by.antibody=TRUE, sel.patterns=NULL){
  # by.antibody = should edges by color-coded by used antibody?
  # sel.nodes = if specified, only a sub-graph is created in which each node names matches to any of these regular expression patterns (see ?grep)
  myGraph <- ftM2graphNEL(assignTab)
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
  myGraphAttrs <- list(node = list(fillcolor=target.color, fontsize="20", fixedsize="FALSE"), edge = list(color = edge.color))
  myNodeAttrs <- list()
  myNodeAttrs$fillcolor <- as.list(antibody.colors)
  myEdges <- edges(myGraph)
  myEdgeAttrs <- list()
  if (by.antibody){
    myEdgeAttrs$color <- list()
    for (this.antibody in my.antibodies)
      for (this.target in myEdges[[this.antibody]])
        myEdgeAttrs$color[[paste(this.antibody,this.target,sep="~")]] <- antibody.colors[this.antibody]
  }#if (by.antibody)
  attr(myGraph,"graphAttrs") <- myGraphAttrs
  attr(myGraph,"nodeAttrs") <- myNodeAttrs
  attr(myGraph,"edgeAttrs") <- myEdgeAttrs
  #class(myGraph) <- c("peakAssignGraph",class(myGraph))
  return(myGraph)
}# assignTable2Graph

myGraphPlot <- function(x,...){
  par(mar=c(0.1,0.1,0.1,0.1), font=2, lwd=2)
  class(x) <- "graphNEL"
  plot(x, attrs=attr(x,"graphAttrs"), nodeAttrs=attr(x,"nodeAttrs"),
       edgeAttrs=attr(x,"edgeAttrs"),...)
}#myGraphPlot

### do testing of functions? might become example code in distant future
if (doFunctionTest) {
  example(findPeaksOnSmoothed)
  assignTabX <- peakList2AssignTable(peaksX,target.name=c("upSymbol","downSymbol"))
  print(assignTabX)
  ### convert table into a graph:
  graphX <- assignTable2Graph(assignTabX)
  myGraphPlot(graphX)
}#if (doFunctionTest)


#############################################################################
### the following is how we used these functions with our TF ChIP-chip data,
### feel free to modify it to suit your own data
#############################################################################

### tell R in which directory the saved binary data resides
dataDir <- "/panfs/panda1/huber/sperling/tfchip/data"

### load list of identified enrichments:
load(file=file.path(dataDir,"peaksV3c.RData"))
## the loaded object allTFPeaks is a list of lists of peaks
tfPeaks <- unlist(allTFPeaks,recursive=FALSE, use.names=FALSE)
if ("modification" %in% names(tfPeaks[[1]])
    tfPeaks <- lapply(tfPeaks, function(p){names(p)[names(p)=="modification"] <- "antibody"; p})
    
assignTF <- peakList2AssignTable(tfPeaks)
assignTF[,1] <- gsub("H1.cdna.2","Dpf3.1", assignTF[,1])
assignTF[,1] <- gsub("H1.cdna.3","Dpf3.2", assignTF[,1])

my.patterns <- unique(c("Nkx2","Myocd", "GATA4","SRF$", "MEF2","bcl2$","Hand","Alk2","Irx4","Tbx","Isl1","Bop","Anf","Dpf3$", paste("^",names(which(table(assignTF[,2])>5)),"$",sep="")))
tf.colors <- brewer.pal(8,"Dark2"); tf.colors[4] <- "#FF00FF"; tf.colors[3] <- "royalblue";
names(tf.colors) <- c("Gata4","Nkx2.5","Mef2","Srf","","Dpf3.1","Dpf3.2","Anti.flag")
graphTF <- assignTable2Graph(assignTF, sel.patterns=my.patterns, antibody.colors=tf.colors)

#pdf(file.path(dataDir,"graphTF.pdf"),width=12,height=10)
myGraphPlot(graphTF)
#dev.off()
