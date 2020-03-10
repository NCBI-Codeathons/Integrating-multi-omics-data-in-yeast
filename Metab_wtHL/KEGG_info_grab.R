# get info for all compounds from KEGG
library(KEGGREST)
# initialize dataframe
cpd.data<-as.data.frame(read.csv('compound_full_list.csv',header=F))
colnames(cpd.data)<-c('ID','Name')

cpd.data$ID<-substring(cpd.data$ID,5,10)
cpd.data$exact.mass<-0
cpd.data$pathway<-'empty'

#cpd.data.list<-split(cpd.data, seq(nrow(cpd.data)))

KeggGrab<-function(cpd.data) {
  print(cpd.data[1])
  temp<-keggGet(cpd.data[1])
  ex.mass<-temp[[1]]$EXACT_MASS
  paw<-paste(names(temp[[1]]$PATHWAY),collapse = ';')
  return(c(ex.mass,paw))
}

cpd.extend<-apply(cpd.data,1,KeggGrab)


for (i in 1:length(cpd.extend)) {
  if(length(cpd.extend[[i]])==1) {
    cpd.extend[[i]][2]<-cpd.extend[[i]][1]
    cpd.extend[[i]][1]<-NA
  }
}

cpd.data$exact.mass<-as.double(sapply(cpd.extend,'[[',1))
cpd.data$pathway<-sapply(cpd.extend,'[[',2)

write.csv(cpd.data,'compound_full_list_anno.csv')
