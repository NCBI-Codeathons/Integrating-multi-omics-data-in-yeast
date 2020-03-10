# convert all compound names to KEGG ID
library(KEGGREST)


cpd.data<-read.csv('mummichog_matched_compound_all.csv')

# matched.Compund is consisted of KEGGID and Biocyc ID (CPD and compound name) 
# and a few compund names
# Biocyc table to convert KEGGID to Biocyc ID and vice versa

#cpd.data$chem.name<-""
#cpd.data$cpd.mass<-NA
cpd.data$KeggID<-"empty"
cpd.data$BiocycID<-'empty'

#cpdname<-as.character(unique(cpd.data$Matched.Compound))

# get keggid
keggidx<-grepl('C[0-9]{5}',cpd.data$Matched.Compound)

keggid<-unique(cpd.data[keggidx,])

cpd.data$KeggID[keggidx]<-c(as.character(cpd.data$Matched.Compound[keggidx]))


# get ID from Biocyc database
cpdidx<-grepl('CPD-[0-9]{1,}',cpd.data$Matched.Compound)

cpd.data$BiocycID[!keggidx]<-c(as.character(cpd.data$Matched.Compound[!keggidx]))

# read in biocyc info
biocyc.cpd.all<-read.csv('Biocyc_all_yeast_compound.csv')

biocyc.cpd.all$ChEBI<-NULL

biocyc.cpd.all$Pathways.of.compound<-gsub("//",";",biocyc.cpd.all$Pathways.of.compound)

# reduce the KEGG column to only KEGG ID
extract.kegg<-function(web.str) {
  web.str<-as.character(web.str)
  if(nchar(web.str)>0) {
    keggreturn<-substr(web.str,nchar(web.str)-9,nchar(web.str)-4)
    return(keggreturn)
  }
}


biocyc.cpd.all$KEGG<-as.character(apply(biocyc.cpd.all[4],1,extract.kegg))

# makes the empty pathway 'empty'
biocyc.cpd.all$Pathways.of.compound[which(nchar(biocyc.cpd.all$Pathways.of.compound)==0)]<-'empty'

# match all available KEGG and Biocyc to make full list of KEGG cpd and biocyc cpd
makeupids<-function(vec.c) {
  kvec<-vec.c[1]
  bvec<-vec.c[2]
  if(kvec=='empty'& is.element(bvec,biocyc.cpd.all$Object.ID)) {
    kvec<-as.character(biocyc.cpd.all$KEGG[which(biocyc.cpd.all$Object.ID==bvec)])
  } else if (bvec=='empty' & is.element(kvec,biocyc.cpd.all$KEGG)) {
    bvec<-as.character(biocyc.cpd.all$Object.ID[which(biocyc.cpd.all$KEGG==kvec)])
  }
  return(list(kvec,bvec))
}


idtemp<-apply(cpd.data[5:6],1,makeupids)

cpd.data$KeggID<-sapply(idtemp,'[[',1)

cpd.data$BiocycID<-sapply(idtemp,'[[',2)

# match cmp directly to KEGG for remaining unannotated ones
# get cpd mass
adduct.w<-read.csv('cpd_adduct_list.csv')

colnames(adduct.w)<-c('form','times','residu')

cal.cpd.mass<-function(cpd.data.s) {
  add.match<-which(adduct.w$form==cpd.data.s[2])
  cpd.mass<-as.numeric(cpd.data.s[1])*adduct.w$times[add.match]+adduct.w$residu[add.match]
  return(cpd.mass)
}

cpd.data$Mass.Diff<-apply(cpd.data[,c(1,3)],1,cal.cpd.mass)
colnames(cpd.data)[4]<-'Cpd.mass'

# read in kegg info
kegg.cpd.all<-read.csv('KEGG_compound_full_list_anno.csv')

chemidx<-(cpd.data$KeggID=='empty' & !cpdidx)

chemidx.n<-c(1:nrow(cpd.data))[chemidx]

cpd.data$chem.name<-'empty'

cpd.data$chem.name[chemidx.n]<-gsub("_|-"," ",cpd.data$Matched.Compound[chemidx.n])


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
    m.dif<-as.numeric(e.mass.vec)-cpd.data$Cpd.mass[i]
    if(min(abs(m.dif),na.rm=T)<2) {
      # has a problem of 2 matching min/cpd, then only get the first one
      cpd.data$KeggID[i]<-substr(names(ktemp[which.min(abs(m.dif))]),5,10)
    }
  }
}


# make col for stats to be before pathway, easier check
cpd.data$p.value<-NA
cpd.data$t.score<-NA

# get kegg pathway and bioccy pathway
cpd.data$keggpath<-'empty'
cpd.data$biocycpath<-'empty'

getkeggpath<-function(kid) {
  if(kid!='empty') {
    kpath<-as.character(kegg.cpd.all$pathway[which(kegg.cpd.all$ID==kid)])
  } else{kpath<-'empty'}
  
  return(kpath)
}

cpd.data$keggpath<-as.character(apply(cpd.data[5],1,getkeggpath))

cpd.data$keggpath[which(cpd.data$keggpath=="character(0)")]<-'empty'

getbiocycpath<-function(bid) {
  if(bid!='empty'){
    bpath<-as.character(biocyc.cpd.all$Pathways.of.compound[which(biocyc.cpd.all$Object.ID==bid)])
  } else {bpath<-'empty'}
  
  return(bpath)
}

cpd.data$biocycpath<-as.character(apply(cpd.data[6],1,getbiocycpath))

cpd.data$biocycpath[which(cpd.data$biocycpath=="character(0)")]<-'empty'

# get peak stats
peak.stats<-read.csv('peak_stats.csv')

cpd.data$Query.Mass<-as.numeric(cpd.data$Query.Mass)

peak.stats$m.z<-as.numeric(peak.stats$m.z)

getpeakstats<-function(pmas){
  matchid<-which(peak.stats$m.z==pmas)
  pval<-peak.stats$p.value[matchid]
  tsc<-peak.stats$t.score[matchid]
  return(list(pval,tsc))
}

stattemp<-apply(cpd.data[1],1,getpeakstats)

cpd.data$p.value<-as.numeric(sapply(stattemp,'[[',1))

cpd.data$t.score<-as.numeric(sapply(stattemp,'[[',2))

write.csv(cpd.data,'wtHL_cpd_path_stats.csv')
