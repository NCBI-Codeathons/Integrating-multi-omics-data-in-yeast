# KEGG integration
library(pathview)
data(gse16873.d)
data(demo.paths)

i<-1

# slow, but output small, original KEGG node labels
pv.out <- pathview(gene.data = gse16873.d[, 1], pathway.id = demo.paths$sel.paths[i],
                   species = "hsa", out.suffix = "gse16873", kegg.native = T)


list.files(pattern="hsa04110", full.names=T)
str(pv.out)
head(pv.out$plot.data.gene)

# fast, output large/layered, node label is offical gene symbol
pv.out <- pathview(gene.data = gse16873.d[, 1], pathway.id = demo.paths$sel.paths[i],
                  species = "hsa", out.suffix = "gse16873.2layer", kegg.native = T,
                  same.layer = F)

# Graphviz, more control over edges and modes
# only change the edge, one layer
pv.out <- pathview(gene.data = gse16873.d[, 1], pathway.id = demo.paths$sel.paths[i],
                   species = "hsa", out.suffix = "gse16873", kegg.native = F,
                   sign.pos = demo.paths$spos[i])

dim(pv.out$plot.data.gene)
head(pv.out$plot.data.gene)

# two layers/pages, 1-main graph; 2-legend
pv.out <- pathview(gene.data = gse16873.d[, 1], pathway.id = demo.paths$sel.paths[i],
                   species = "hsa", out.suffix = "gse16873.2layer", kegg.native = F,
                   sign.pos = demo.paths$spos[i], same.layer = F)

# split the node groups into individual detached nodes
pv.out <- pathview(gene.data = gse16873.d[, 1], pathway.id = demo.paths$sel.paths[i],
                   species = "hsa", out.suffix = "gse16873.split", kegg.native = F,
                   sign.pos = demo.paths$spos[i], split.group = T)
dim(pv.out$plot.data.gene)

# expand the multiple-gene nodes into individual genes
pv.out <- pathview(gene.data = gse16873.d[, 1], pathway.id = demo.paths$sel.paths[i],
                   species = "hsa", out.suffix = "gse16873.split.expanded", kegg.native = F,
                   sign.pos = demo.paths$spos[i], split.group = T, expand.node = T)
dim(pv.out$plot.data.gene)
head(pv.out$plot.data.gene)
# split nodes or expanded genes may inherit the edges from the unsplit group or unexpanded nodes
# get a gene/protein-gene/protein interaction network
# better view the network characteristics (modularity etc) and gene-wise (instead of node-wise) data.
# instead of ndoe-wise data


###############
# data integration
###############

# simulate some compound data
sim.cpd.data=sim.mol.data(mol.type="cpd", nmol=3000)
data(cpd.simtypes)

i <- 3
print(demo.paths$sel.paths[i])

# generate KEGG view of a pathway with gene and compound
pv.out <- pathview(gene.data = gse16873.d[, 1], cpd.data = sim.cpd.data,
                   pathway.id = demo.paths$sel.paths[i], species = "hsa", out.suffix = "gse16873.cpd",
                   keys.align = "y", kegg.native = T, key.pos = demo.paths$kpos1[i])

str(pv.out)

head(pv.out$plot.data.cpd)
# compound without data is filled with white

#plot Graphviz, better hierarchical structure
#reaction is parsed into relatioship between gene and compound
pv.out <- pathview(gene.data = gse16873.d[, 1], cpd.data = sim.cpd.data,
                   pathway.id = demo.paths$sel.paths[i], species = "hsa", out.suffix = "gse16873.cpd",
                   keys.align = "y", kegg.native = F, key.pos = demo.paths$kpos2[i],
                   sign.pos = demo.paths$spos[i], cpd.lab.offset = demo.paths$offs[i])

#######
#integrate multiple samples and plot
#######

# simulate compound data
set.seed(10)
sim.cpd.data2 = matrix(sample(sim.cpd.data, 18000,replace = T), ncol = 6)
rownames(sim.cpd.data2) = names(sim.cpd.data)
colnames(sim.cpd.data2) = paste("exp", 1:6, sep = "")
head(sim.cpd.data2, 3)

# gene.data has 3 samples and cpd.data has 2.
# gene and cpd nodes are sliced into multiple pieces according to samples
# as samples size is different between gene and cpd, we can match the paired samples
# paired samples--gene and cpd comes from the same sample

# KEGG view without data matching
pv.out <- pathview(gene.data = gse16873.d[, 1:3],
                   cpd.data = sim.cpd.data2[, 1:2], pathway.id = demo.paths$sel.paths[i],
                   species = "hsa", out.suffix = "gse16873.cpd.3-2s", keys.align = "y",
                   kegg.native = T, match.data = F, multi.state = T, same.layer = T)
head(pv.out$plot.data.cpd)

# KEGG with data matching
pv.out <- pathview(gene.data = gse16873.d[, 1:3],
                   cpd.data = sim.cpd.data2[, 1:2], pathway.id = demo.paths$sel.paths[i],
                   species = "hsa", out.suffix = "gse16873.cpd.3-2s.match",
                   keys.align = "y", kegg.native = T, match.data = T, multi.state = T,
                   same.layer = T)

#graphviz view without data matching
pv.out <- pathview(gene.data = gse16873.d[, 1:3],
                   cpd.data = sim.cpd.data2[, 1:2], pathway.id = demo.paths$sel.paths[i],
                   species = "hsa", out.suffix = "gse16873.cpd.3-2s", keys.align = "y",
                   kegg.native = F, match.data = F, multi.state = T, same.layer = T)

# plot sample separately on each graph, gene and cpd has to match to be on the same page
# for sample with only one type of data, the other type will be all white
pv.out <- pathview(gene.data = gse16873.d[, 1:3],
                   cpd.data = sim.cpd.data2[, 1:2], pathway.id = demo.paths$sel.paths[i],
                   species = "hsa", out.suffix = "gse16873.cpd.3-2s", keys.align = "y",
                   kegg.native = T, match.data = F, multi.state = F, same.layer = T)
# KEGG on two layers to speed up, then the label is not KEGG 
pv.out <- pathview(gene.data = gse16873.d[, 1:3],
                   cpd.data = sim.cpd.data2[, 1:2], pathway.id = demo.paths$sel.paths[i],
                   species = "hsa", out.suffix = "gse16873.cpd.3-2s.2layer",
                   keys.align = "y", kegg.native = T, match.data = F, multi.state = T,
                   same.layer = F)

#######
# discrete data
#######
# with genes or cpds that are significant (1 or 0) or shorter list 
require(org.Hs.eg.db)
gse16873.t <- apply(gse16873.d, 1, function(x) t.test(x,alternative = "two.sided")$p.value)
sel.genes <- names(gse16873.t)[gse16873.t < 0.1]
sel.cpds <- names(sim.cpd.data)[abs(sim.cpd.data) > 0.5]

# plot both gene and cpd as discrete numbers
# only select 0<gene<5 and 0<cpd<2
# have 5 bins for gene and 2 bins for cpd
pv.out <- pathview(gene.data = sel.genes, cpd.data = sel.cpds,
                     pathway.id = demo.paths$sel.paths[i], species = "hsa", out.suffix = "sel.genes.sel.cpd",
                     keys.align = "y", kegg.native = T, key.pos = demo.paths$kpos1[i],
                     limit = list(gene = 5, cpd = 2), bins = list(gene = 5, cpd = 2),
                     na.col = "gray", discrete = list(gene = T, cpd = T))

# plot only gene as discrete but cpd as continuous
# select 0<gene <5 but -1<cpd<1
# have 5 bins for gene and 10 binds for cpd
pv.out <- pathview(gene.data = sel.genes, cpd.data = sim.cpd.data,
                     pathway.id = demo.paths$sel.paths[i], species = "hsa", out.suffix = "sel.genes.cpd",
                     keys.align = "y", kegg.native = T, key.pos = demo.paths$kpos1[i],
                     limit = list(gene = 5, cpd = 1), bins = list(gene = 5, cpd = 10),
                     na.col = "gray", discrete = list(gene = T, cpd = F))

########
# ID mapping
########

# map external ID types to standard KEGG IDs -- automatically
# specify the external ID types using gene.idtype and cpd.idtype arguments
cpd.cas <- sim.mol.data(mol.type = "cpd", id.type = cpd.simtypes[2],nmol = 10000)
gene.ensprot <- sim.mol.data(mol.type = "gene", id.type = gene.idtype.list[4],nmol = 50000)

pv.out <- pathview(gene.data = gene.ensprot, cpd.data = cpd.cas,
                   gene.idtype = gene.idtype.list[4], cpd.idtype = cpd.simtypes[2],
                   pathway.id = demo.paths$sel.paths[i], species = "hsa", same.layer = T,
                   out.suffix = "gene.ensprot.cpd.cas", keys.align = "y", kegg.native = T,
                   key.pos = demo.paths$kpos2[i], sign.pos = demo.paths$spos[i],
                   limit = list(gene = 3, cpd = 3), bins = list(gene = 6, cpd = 6))

# manually create external ID and KEGG ID matrix
id.map.cas <- cpdidmap(in.ids = names(cpd.cas), in.type = cpd.simtypes[2],
                       out.type = "KEGG COMPOUND accession")
cpd.kc <- mol.sum(mol.data = cpd.cas, id.map = id.map.cas)
id.map.ensprot <- id2eg(ids = names(gene.ensprot),category = gene.idtype.list[4], org = "Hs")
gene.entrez <- mol.sum(mol.data = gene.ensprot, id.map = id.map.ensprot)
pv.out <- pathview(gene.data = gene.entrez, cpd.data = cpd.kc,
                   pathway.id = demo.paths$sel.paths[i], species = "hsa", same.layer = T,
                   out.suffix = "gene.entrez.cpd.kc", keys.align = "y", kegg.native = T,
                   key.pos = demo.paths$kpos2[i], sign.pos = demo.paths$spos[i],
                   limit = list(gene = 3, cpd = 3), bins = list(gene = 6, cpd = 6))
#############
# different species
############

data(korg)
head(korg)

#number of species which use Entrez Gene as the default ID
sum(korg[,"entrez.gnodes"]=="1")

#number of species which use other ID types or none as the default ID
sum(korg[,"entrez.gnodes"]=="0")

#species which do not have Entrez Gene annotation at all
na.idx=is.na(korg[,"ncbi.geneid"])
head(korg[na.idx,])

# species you can work with other ID types
data(bods)
bods
data(gene.idtype.list)
gene.idtype.list

# work with other ID types
sce.dat.entrez <- sim.mol.data(mol.type="gene",id.type="entrez",species="sce",nmol=3000)

pv.out <- pathview(gene.data = sce.dat.entrez, gene.idtype="entrez",
                   pathway.id = "00640", species = "sce", out.suffix = "sce.entrez",
                   kegg.native = T, same.layer=T)

egid.sce=eg2id(names(sce.dat.entrez), category="GENENAME", pkg="org.Sc.sgd.db")
sce.dat.gn <- sce.dat.entrez
names(sce.dat.gn) <- egid.sce[,2]
head(sce.dat.gn)

pv.out <- pathview(gene.data = sce.dat.gn, gene.idtype="genename",
                   pathway.id = "00640", species = "sce", out.suffix = "sce.genename.2layer",
                   kegg.native = T, same.layer=F)

# analysis of unannotated species with gene mapped to orthologs
# species='ko'
ko.data=sim.mol.data(mol.type="gene.ko", nmol=5000)
pv.out <- pathview(gene.data = ko.data, pathway.id = "04112",
                   species = "ko", out.suffix = "ko.data", kegg.native = T)


