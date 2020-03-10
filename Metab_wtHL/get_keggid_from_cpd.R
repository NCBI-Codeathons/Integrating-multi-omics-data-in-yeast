# get chemical compound
chemidx<-!keggidx & !cpdidx

chemidx.n<-c(1:nrow(cpd.data))[chemidx]

cpd.data$chem.name[chemidx.n]<-gsub("_|-"," ",cpd.data$Matched.Compound[chemidx.n])


#write.csv(cpd.data,file='cpd_name_forKEGGID.csv')

adduct.list<-unique(cpd.data$Matched.Form)

write.csv(adduct.list,file='cpd_adduct_list_raw.csv')

# get chem mass

adduct.w<-read.csv('cpd_adduct_list.csv')

colnames(adduct.w)<-c('form','times','residu')

cal.cpd.mass<-function(cpd.data.s) {
  add.match<-which(adduct.w$form==cpd.data.s[2])
  cpd.mass<-as.numeric(cpd.data.s[1])*adduct.w$times[add.match]+adduct.w$residu[add.match]
  return(cpd.mass)
}

cpd.data$cpd.mass<-apply(cpd.data[,c(1,3)],1,cal.cpd.mass)

# readin all kegg info

kegg.cpd.all<-read.csv('compound_full_list_anno.csv')

# get kegg id for compound
for (i in chemidx.n) {
  ktemp<-keggFind('compound',cpd.data$chem.name[i])
  print(i)
  if (length(ktemp)==1) {
    cpd.data$KeggID[i]<-substr(names(ktemp),5,10)
  } else if (length(ktemp)>1) {
    e.mass.vec<-vector(mode='numeric',length=length(ktemp))
    for (n in 1:length(ktemp)) {
      print(paste(i,n,sep='-'))
      ktemp.n<-substr(names(ktemp[n]),5,10)
      e.mass.vec[n]<-max(kegg.cpd.all$exact.mass[which(kegg.cpd.all$ID==ktemp.n)],0,na.rm=T)
    }
    m.dif<-as.numeric(e.mass.vec)-cpd.data$cpd.mass[i]
    if(min(abs(m.dif),na.rm=T)<2) {
      # has a problem of 2 matching min/cpd, then only get the first one
      cpd.data$KeggID[i]<-substr(names(ktemp[which.min(abs(m.dif))]),5,10)
    }
  }
}
