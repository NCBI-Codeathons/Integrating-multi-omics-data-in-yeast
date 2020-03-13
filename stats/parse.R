parseReaction <- function(reaction) {
  attrs <- xmlAttrs(reaction)
  
  ## required: name,type
  name <- attrs[["name"]]
  type <- attrs[["type"]]
  
  children <- xmlChildren(reaction)
  
  ## more than one substrate/product possible
  childrenNames <- names(children)
  substrateIndices <- grep("^substrate$", childrenNames)
  productIndices <- grep("^product$", childrenNames)
  substrateName <- substrateAltName <- vector("character", length(substrateIndices))
  productName <- productAltName <- vector("character", length(productIndices))  
  
  for (i in seq(along=substrateIndices)) {
    ind <- substrateIndices[i]
    substrate <- children[[ind]]
    substrateName[i] <- xmlAttrs(substrate)[["name"]]
    substrateAltName[i] <- as.character(NA)
    
    substrateChildren <- .xmlChildrenWarningFree(substrate)
    if (!is.null(substrateChildren)) {
      substrateAlt <- substrateChildren$alt
      substrateAltName[i] <- xmlAttrs(substrateAlt)[["name"]]
    }
    
  }
  
  for(i in seq(along=productIndices)) {
    ind <- productIndices[i]
    product <- children[[ind]]
    productName[i] <- xmlAttrs(product)[["name"]]
    productChildren <- .xmlChildrenWarningFree(product)
    productAltName[i] <- as.character(NA)
    if(!is.null(productChildren)) {
      productAlt <- productChildren$alt
      productAltName[i] <- xmlAttrs(productAlt)[["name"]]
    }
  }
  
  new("KEGGReaction",
      name = name,
      type = type,
      substrateName = substrateName,
      substrateAltName = substrateAltName,
      productName = productName,
      productAltName = productAltName)
}
parseEntry <- function(entry) {
  attrs <- xmlAttrs(entry)
  
  ## required: id, name,type
  entryID <- attrs[["id"]]
  name <- unname(unlist(strsplit(attrs["name"]," ")))
  type <- attrs[["type"]]
  
  ## implied: link, reaction, map
  link <- getNamedElement(attrs,"link")
  reaction <- getNamedElement(attrs, "reaction")
  map <- getNamedElement(attrs, "map")
  
  ## graphics
  graphics <- xmlChildren(entry)$graphics
  g <- parseGraphics(graphics)
  
  ## types: ortholog, enzyme, gene, group, compound and map
  if(type != "group") {
    newNode <- new("KEGGNode",
                   entryID=entryID,
                   name=name,
                   type=type,
                   link=link,
                   reaction=reaction,
                   map=map,
                   graphics=g)
  } else if(type=="group") {
    children <- xmlChildren(entry)
    children <- children[names(children) == "component"]
    if(length(children)==0) {
      component <- as.character(NA)
    } else {
      component <- sapply(children, function(x) {
        if(xmlName(x) == "component") {
          return(xmlAttrs(x)["id"])
        } else {
          return(as.character(NA))
        }     
      })
    }
    component <- unname(unlist(component))
    newNode <- new("KEGGGroup",
                   component=component,
                   entryID=entryID,
                   name=name,
                   type=type,
                   link=link,
                   reaction=reaction,
                   map=map,
                   graphics=g
    )
  }
  return(newNode)
}
parseRelation <- function(relation) {
  attrs <- xmlAttrs(relation)
  
  ## required: entry1, entry2, type
  entry1 <- attrs[["entry1"]]
  entry2 <- attrs[["entry2"]]
  type <- attrs[["type"]]
  
  subtypeNodes <- xmlChildren(relation)
  subtypes <- sapply(subtypeNodes, parseSubType)
  newEdge <- new("KEGGEdge",
                 entry1ID=entry1,
                 entry2ID=entry2,
                 type=type,
                 subtype=subtypes
  )                     
  return(newEdge)
}
parsePathwayInfo <- function(root) {
  attrs <- xmlAttrs(root)
  ## required: name, org, number
  name <- attrs[["name"]]
  org <- attrs[["org"]]
  number <- attrs[["number"]]
  ## implied: title, image, link
  title <- getNamedElement(attrs, "title")
  image <- getNamedElement(attrs, "image")
  link <- getNamedElement(attrs, "link")
  
  return(new("KEGGPathwayInfo",
             name=name,
             org=org,
             number=number,
             title=title,
             image=image,
             link=link))
}
parseGraphics <- function(graphics) {
  if(is.null(graphics))
    return(new("KEGGGraphics"))
  attrs <- xmlAttrs(graphics)
  g <- new("KEGGGraphics",
           name=getNamedElement(attrs,"name"),
           x=as.integer(getNamedElement(attrs,"x")),
           y=as.integer(getNamedElement(attrs,"y")),
           type=getNamedElement(attrs,"type"),
           width=as.integer(getNamedElement(attrs, "width")),
           height=as.integer(getNamedElement(attrs,"height")),
           fgcolor=getNamedElement(attrs, "fgcolor"),
           bgcolor=getNamedElement(attrs, "bgcolor")
  )
  return(g)
  
}
parseSubType <- function(subtype) {
  attrs <- xmlAttrs(subtype)
  name <- attrs[["name"]]
  value <- attrs[["value"]]
  return(new("KEGGEdgeSubType",name=name, value=value))
}
getNamedElement <- function(vector, name) {
  if (name %in% names(vector))
    return(vector[[name]])
  else
    return(as.character(NA))
}
.xmlChildrenWarningFree <- function(xmlNode) {
  if(is.null(xmlNode$children))
    return(NULL)
  return(XML::xmlChildren(xmlNode))
}
library("KEGGREST") 
total_path.l<-unique(keggLink("pathway","sce"))
path_id<-sapply(total_path.l,function(i) strsplit(i,split=":")[[1]][2]) 
sapply(path_id,function(id){
  fn<-paste(id,".xml",sep="",collapse = "")
  print(fn)
  doc <- xmlTreeParse(fn, getDTD=FALSE,error=xmlErrorCumulator(immediate=FALSE))
  r <- xmlRoot(doc)
  
  ## possible elements: entry, relation and reaction
  childnames <- sapply(xmlChildren(r), xmlName)
  isEntry <- childnames == "entry"
  isRelation <- childnames == "relation"
  isReaction <- childnames == "reaction"
  
  ## parse them
  kegg.pathwayinfo <- parsePathwayInfo(r)
  kegg.nodes <- sapply(r[isEntry], parseEntry)
  kegg.edges <- sapply(r[isRelation], parseRelation)
  kegg.reactions <- sapply(r[isReaction], parseReaction)
  names(kegg.nodes) <- sapply(kegg.nodes, getEntryID)
  
  
  gene_rn.m<-Reduce(rbind,sapply(kegg.nodes,function(i){
    if(length(i@name)>=1){
      if(grepl("sce:",i@name[1])){
        cbind(i@name,rep(i@reaction,length(i@name)))
      }
    }
  }))
  meta_rn.m<-Reduce(rbind,lapply(kegg.reactions,function(i){
    if(length(i@substrateName)>=1){
      re1<-cbind(i@substrateName,rep(i@name,length(i@substrateName)))
    }
    if(length(i@productName)>=1){
      re2<-cbind(i@productName,rep(i@name,length(i@productName)))
      re<-rbind(re1,re2)
    }
  }))
  if(!is.null(dim(meta_rn.m))){
    gene_meta.m<-Reduce(rbind,sapply(unique(gene_rn.m[,2]),function(r){
      g.idx<-which(gene_rn.m[,2]==r)
      m.idx<-which(meta_rn.m[,2]==r)
      if(length(c(g.idx,m.idx))>max(length(m.idx),length(g.idx))){
        re<-sapply(gene_rn.m[g.idx,1], function(a) lapply(meta_rn.m[m.idx,1], function (b) c(as.character(a), as.character(b))) )
        matrix(unlist(re),ncol=2,byrow = T)
      }
  }))
    out.fn<-paste(id,"_gm_edge.csv",sep = "",collapse = "")
    write.csv(gene_meta.m,out.fn,row.names=F)
  }
}) 
file="sce00040.xml"
doc <- xmlTreeParse(file, getDTD=FALSE,error=xmlErrorCumulator(immediate=FALSE))
r <- xmlRoot(doc)

## possible elements: entry, relation and reaction
childnames <- sapply(xmlChildren(r), xmlName)
isEntry <- childnames == "entry"
isRelation <- childnames == "relation"
isReaction <- childnames == "reaction"

## parse them
kegg.pathwayinfo <- parsePathwayInfo(r)
kegg.nodes <- sapply(r[isEntry], parseEntry)
kegg.edges <- sapply(r[isRelation], parseRelation)
kegg.reactions <- sapply(r[isReaction], parseReaction)
names(kegg.nodes) <- sapply(kegg.nodes, getEntryID)


gene_rn.m<-Reduce(rbind,sapply(kegg.nodes,function(i){
  if(length(i@name)>=1){
    if(grepl("sce:",i@name[1])){
      cbind(i@name,rep(i@reaction,length(i@name)))
    }
  }
}))
meta_rn.m<-Reduce(rbind,sapply(kegg.reactions,function(i){
  if(length(i@substrateName)>=1){
    re1<-cbind(i@substrateName,rep(i@name,length(i@substrateName)))
  }
  if(length(i@productName)>=1){
    re2<-cbind(i@productName,rep(i@name,length(i@productName)))
    re<-rbind(re1,re2)
  }
}))
Reduce(rbind,sapply(unique(gene_rn.m[,2]),function(r){
  g.idx<-which(gene_rn.m[,2]==r)
  m.idx<-which(meta_rn.m[,2]==r)
  re<-sapply(gene_rn.m[g.idx,1], function(a) lapply(meta_rn.m[m.idx,1], function (b) c(as.character(a), as.character(b))) )
  matrix(unlist(re),ncol=2,byrow = T)
}))
