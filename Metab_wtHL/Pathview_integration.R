# KEGG integration
library(pathview)

# load RNA-seq data
head(seqdata.norm)

# load compound data

#rm(list=setdiff(ls(), c('seqdata.norm','gene.idtype.bods','korg','cpd.simtypes','gene.idtype.list')))

#detach("package:pathview", unload=TRUE)
#library(pathview)


# plot gene data for a specific pathway
pv.out <- pathview(gene.data = seqdata.norm[, 1], pathway.id = '00730', 
                   gene.idtype='ORF',species = "sce", out.suffix = "thiamin", kegg.native = T)

pv.out <- pathview(gene.data = seqdata.norm[, 1], pathway.id = '00730', 
                   gene.idtype='ORF',species = "sce", out.suffix = "thiamin_graph", kegg.native = F,sign.pos = 'bottomright')


