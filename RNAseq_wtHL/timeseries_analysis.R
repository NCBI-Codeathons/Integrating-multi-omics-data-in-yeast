############
# functions for later
############
cols <- function(lowi = "green", highi = "red", ncolors = 20) {
  low <- col2rgb(lowi)/255
  high <- col2rgb("black")/255
  col1 <- rgb(seq(low[1], high[1], len = ncolors), seq(low[2], 
                                                       high[2], len = ncolors), seq(low[3], high[3], len = ncolors))
  low <- col2rgb("black")/255
  high <- col2rgb(highi)/255
  col2 <- rgb(seq(low[1], high[1], len = ncolors), seq(low[2], 
                                                       high[2], len = ncolors), seq(low[3], high[3], len = ncolors))
  col<-c(col1[1:(ncolors-1)],col2)
  return(col)
}

hclust2 <- function(x, method="average", ...)
  hclust(x, method=method, ...)
dist2 <- function(x, ...)
  as.dist(1-cor(t(x), method="pearson"))

##############
##############
dataMatrixFile<- "counts_all.txt"
metaDataFile<- "meta.txt"

# Read dataMatrixFile
seqdata.raw<-read.table(dataMatrixFile,sep="\t",header=T,row.names=1,fill=T)
summary(seqdata)
seqdata<-seqdata.raw
# Strip Annotation and reassign the x matrix

ann<- seqdata[,1:5]
seqdata<- seqdata[,-1*1:5]

# Read metaDataFile

expmeta<-read.table(metaDataFile,sep="\t",header=T,row.names=1)

seqdata<- seqdata[,rownames(expmeta)]
gene.sample.names<- dimnames(seqdata)

seqdata<- apply(seqdata,2,function(x){as.numeric(as.vector(x))})
sample.sum<- round( colSums(seqdata) / 1e6, 1 )

dimnames(seqdata)<- gene.sample.names
colnames(seqdata)<- gsub("_sorted.bam","",colnames(seqdata))
rownames(expmeta)<- gsub("_sorted.bam","",rownames(expmeta))


library(DESeq2)
# Normalize--Variance-stabilizing transformation
# The aim behind the choice of a variance-stabilizing transformation is to find a simple function ƒ to apply to values x 
# in a data set to create new values y = ƒ(x) such that the variability of the values y is not related to their mean value.

# The point of these transforms is to reduce (ideally eliminate) dependence of the variance on the mean.
# particularly the high variance of the logarithm of count data when the mean is low

# Many common statistical methods for exploratory analysis of multidimensional data, 
# especially methods for clustering and ordination (e.g., principal-component analysis and the like), 
# work best for (at least approximately) homoskedastic data; 
# this means that the variance of an observable quantity (i.e., here, the expression strength of a gene) does not depend on the mean.

# You should prefer to use these when you are doing downstream analysis on your count data that doesn't involve testing for 
# differential expression using the statistical methods developed for count data. 
# These scenarios include doing things like clustering, 
# or PCA over your expression data or using the data as input to another machine learning algorithm
seqdata.norm <- vst(seqdata)

boxplot(seqdata.norm+1)
# TAF10, YDR167W, is housekeeping gene
# YNL219C, ALG9
# ALG9,TAF10, TFC1, UBC6 and to a lesser extent KRE11
# should be stable across samples
barplot(seqdata[rownames(ann)=="YNL219C",],las=2,ylab="YNL219C")
barplot(seqdata.norm[rownames(ann)=="YNL219C",],las=2,ylab="YNL219C")
barplot(seqdata[rownames(ann)=="YDR167W",],las=2,ylab="YDR167W")
barplot(seqdata.norm[rownames(ann)=="YDR167W",],las=2,ylab="YDR167W")
#boxplot(normFactor~expmeta$timepoint)
boxplot(sample.sum~expmeta$timepoint)

# and row center for later visualization
# centering is done by subtracting the column means (omitting NAs) of x from their corresponding columns
# scale is TRUE then scaling is done by dividing the (centered) columns of x by their standard deviations if center is TRUE,
# and the root mean square otherwise.
seqdata.norm.genecenter<- t(scale(t(seqdata.norm),scale=F))

chrs<- substr(ann[,1],0,5)
chrs<- gsub(";.*","",chrs)

# filter low counts
filt<- rowSums(seqdata)>10
seqdata<- seqdata[filt,]
seqdata.norm<- seqdata.norm[filt,]
seqdata.norm.genecenter<- seqdata.norm.genecenter[filt,]
ann<- ann[filt,]
chrs<- chrs[filt]
seqdata.raw<-seqdata.raw[filt,]

#############
# data visualization with clustering
############
sds<- apply(seqdata.norm,1,sd,na.rm=T)
hist(sds)
# sds as function of chr
boxplot(sds~chrs)
table(sds>2.75)

# scale mean but not variation
heatmap(seqdata.norm.genecenter[sds>0.4,],col=cols(),hclustfun=hclust2,distfun=dist2)
heatmap(seqdata.norm.genecenter[sds>0.4,],col=cols(),hclustfun=hclust2,distfun=dist2,scale='none')

label<- matrix(c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#a65628','#f781bf','#999999')[as.numeric(factor(expmeta$timepoint))])
heatmap(seqdata.norm.genecenter[sds>3,],col=cols(),hclustfun=hclust2,distfun=dist2,ColSideColors=label)
heatmap(seqdata.norm.genecenter[sds>2.75,],col=cols(),hclustfun=hclust2,distfun=dist2,scale='none')

library(gplots)
temp<-seqdata.norm.genecenter[sds>2.75,]
heatmap.2(temp,dendrogram = 'col',col=cols(),hclustfun=hclust2,distfun=dist2,scale='none',
          trace='none',density.info = 'none',keysize=1)

pdf("sds275.pdf",width=10,height=8)
temp<-seqdata.norm.genecenter[sds>2.75,]
heatmap.2(temp,dendrogram = 'col',col=cols(),hclustfun=hclust2,distfun=dist2,scale='none',
          trace='none',density.info = 'none',keysize=1)
dev.off()




table(sds>2)
gfilt.all<-sds>2
genelist.all<-seqdata[gfilt.all,]
write.table(genelist.all, file='genelist_sds2.txt', sep='\t',row.names=T,col.names=T,quote=F)


###
# compare wt00 and wt02
###

# convert to RPKM

totalR<-colSums(seqdata)
totalRM<-totalR/1000000

seqdata_RPM<-t(t(seqdata)/totalRM)

seqdata_RPKM<-seqdata_RPM/(seqdata.raw$Length/1000)

logf02<-log2(seqdata_RPKM[,2]/seqdata_RPKM[,1])

table(abs(logf02)>3)


lfilt<-abs(logf02)>3
lfilt[is.na(lfilt)]<-FALSE
temp<-seqdata.norm.genecenter[lfilt,]
heatmap.2(tx,dendrogram = 'col',col=cols(),hclustfun=hclust2,distfun=dist2,scale='none',
          trace='none',density.info = 'none',keysize=1)

write.table(temp, file='wt00vs02.txt', sep='\t',row.names=T,col.names=T,quote=F)


#############
# data visualization with pca
############
pca<- prcomp(t(seqdata.norm.genecenter))
plot(pca$x[,1],pca$x[,2],col=label)
tp<- matrix(c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#a65628','#f781bf','#999999')[as.numeric(factor(expmeta$timepoint))])
plot(pca$x[,1],pca$x[,2],col=tp,pch=19,xlab="PC1",ylab="PC2")
# proportion of each component contribution
pairs(pca$x[,1:4],col=tp,pch=19)
barplot(pca$sdev)
# proportion of the variance that each eigenvector represents
# can be calculated by dividing the eigenvalue corresponding to that eigenvector by the sum of all eigenvalues
signif(pca$sdev/sum(pca$sdev)*100,3)
plot(pca$x[,1],log10(sample.sum),xlab="PC1",col=tp,pch=19,ylab="Log10 total reads")

# pca of normalized and gene centered data
pca<- prcomp(seqdata.norm.genecenter)
plot(pca$x[,1],pca$x[,2],col=tp,pch=19,xlab="PC1",ylab="PC2")
# plto pairs of pca components
pairs(pca$x[,1:4],col=tp,pch=19)

# proportion of each component contribution
barplot(pca$sdev,names=paste("PC",1:dim(pca$x)[[2]],sep=""))
signif(pca$sdev/sum(pca$sdev)*100,3)




# write results
write.table(results, file="deseq.results.txt", sep="\t",col.names=NA)
write.table(results[gfilt,], file="deseq.results_sig.txt", sep="\t",col.names=NA)
