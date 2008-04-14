### use biomaRt to construct an annotation data frame for all genomic features
## for a species, as currently annotated in the Ensembl data base

library("biomaRt")
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
allChr <- c(1:22,"X","Y")

gene.ids <- unique(unlist(lapply(as.list(allChr), function(this.chr) getBM(attributes="ensembl_gene_id", filters="chromosome_name", values=this.chr, mart=ensembl)[,1]), use.names=FALSE))

# what information to get for each transcript:
sel.attributes=c("ensembl_gene_id", "ensembl_transcript_id", "hgnc_symbol", "chromosome_name", "strand", "start_position","end_position", "transcript_start", "transcript_end", "description")
# if 'hgnc_symbol' is not a defined attribute for gene symbols in your data species, try to find an equivalent, using commands such as this one:
# grep("symbol",listAttributes(ensembl)[,1], value=TRUE)

# retreive information:
gff <- getBM(attributes=sel.attributes, filters="ensembl_gene_id", value=gene.ids, mart=ensembl)

## replace attribute names by standardized names
gff$gene <- gff$"ensembl_gene_id"
gff$name <- gff$ensembl_transcript_id
gff$chr <- gff$chromosome_name
gff$symbol <- gff$"hgnc_symbol"
gff$start <- gff$transcript_start
gff$end <- gff$transcript_end
gff$"gene_start" <- gff$"start_position"
gff$"gene_end" <- gff$"end_position"
gff$feature <- rep("transcript",nrow(gff))
gff$chromosome_name <- gff$markersymbol <- NULL
gff$transcript_start <- gff$transcript_end <- NULL
gff$start_position <- gff$end_position <- NULL
gff$ensembl_transcript_id <- NULL
gff$ensembl_gene_id <- NULL

  
### some transcripts names my occurr in multiples, usually this is because an Ensembl transcript can have more than one Symbol defined for it:
if (any(duplicated(gff$name))){
  dupl.trans <- unique(gff$name[duplicated(gff$name)])
  G <- lapply(as.list(dupl.trans), function(this.trans){
    this.gff <- subset(gff,name == this.trans)
    if (nrow(unique(this.gff[,c("gene","name","chr","start","end","description")]))>1) return(this.gff[1,,drop=FALSE])
    non.zero.gff <- subset(this.gff, nchar(symbol)>0)
    this.other.sym <- NULL
    if (nrow(non.zero.gff)> 0){
      this.new.sym <- non.zero.gff$symbol[which.min(nchar(non.zero.gff$symbol))]
      if (nrow(non.zero.gff)>1)
        this.other.sym <- paste("Synonyms",paste(non.zero.gff$symbol[-which.min(nchar(non.zero.gff$symbol))],collapse=","),sep=":")
    } else { this.new.sym <- "" }
    this.gff$symbol[1] <- this.new.sym
    if (!is.null(this.other.sym))
      this.gff$description[1] <- paste(this.gff$description[1],this.other.sym,sep=";")
    return(this.gff[1,,drop=FALSE])
  })
  gff <- rbind(gff[-which(gff$name %in% dupl.trans),],do.call("rbind",G))
  rm(G, dupl.trans);gc()
}#if (any(duplicated(gff$name)))


# reorder by chromosome and start position
gff <- gff[order(gff$chr, gff$gene_start, gff$start),c("name","gene","chr","strand","gene_start","gene_end","start","end","symbol", "description","feature")]
rownames(gff) <- NULL
  
# save(gff, file="mygff.RData")
