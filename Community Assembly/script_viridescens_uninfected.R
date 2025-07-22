library(vegan);library(parallel);library(microbiome);library(picante);library(ape)

# COMMUNITY ASSEMBLY ANALYSIS

# Part 1: Calculate bNTI

#Load files

#1-Load ASV table as matrix
viridescens.otu.table <- read.csv("newts_otu_table_rarefy.csv", header = T, check.names = F, row.names = 1)
viridescens.otu.matrix <- as.matrix(viridescens.otu.table)

samp.group <- read.table("newts_metadata.csv", head = T, sep = ",")
samp.meta <- read.table("newts_metadata.csv", head = T, sep = ",")

#Load phylogenetic tree
library(picante)
viridescens.tree <- ape::read.tree(file="newts_final_root.nwk")
viridescens.tree <- ape::multi2di(viridescens.tree)

# Calculate distance
dis <- cophenetic.phylo(viridescens.tree)# calculate phylogenetic distance matrix
spname <- rownames(dis)# names of species with phylogenetic distance
colnames(dis) <- rownames(dis)

# Match the data
comm <- viridescens.otu.matrix[match(spname,rownames(viridescens.otu.matrix)), 1:ncol(viridescens.otu.matrix)]# remove the species without phylogenetic information.
rownames(comm) <- spname
comm <- comm[,match(samp.meta[,1],colnames(comm))]# sort community data by sample grouping file. 
comm <- t(comm) # need transposition

source("../bNTI.R")# load function of betaNTI, for Hydra change this to full code

# Observed betaMNTD across all samples
bMNTD<-comdistnt(comm, dis, abundance.weighted = TRUE, exclude.conspecifics = FALSE) # calculate observed betaMNTD. pay attention to weight by abundance or not.
write.csv(as.matrix(bMNTD),file="bMNTD.csv")# output observed betaMNTD

# Calculate betaNTI across all samples
bNTI.all=bNTI(comm, dis, samp.group=NA, weighted=TRUE,grouping=FALSE,rand=1000,output.bMNTD=FALSE)
write.csv(bNTI.all,file="bNTIall.csv")# output

# Calculate betaMNTD and betaNTI within group
bNTI.group=bNTI(comm, dis, samp.group=samp.group, weighted=TRUE,grouping=TRUE,rand=1000,output.bMNTD=TRUE)

write.csv(bNTI.group$betaNTI,file="bNTI.group_viridescens.csv")
write.csv(bNTI.group$betaMNTD,file="bMNTD.group_viridescens.csv")

# Part 2: 
com <- read.csv("newts_uninfected_otu_table_rare.csv", header = T, check.names = F, row.names = 1)
comm <- com[,colSums(com)>0] # remove undetected species
comm <- as.data.frame(t(comm)) # Transpose dataframe

# Calculate Bray Curtis
BC <- vegdist(comm,method="bray")
write.csv(as.matrix(BC),"BC.csv")

# Calculate RC
source("../RC.p.R")
rc<-RC.p(comm,method="bray",rand=100,portion=FALSE,
         nworker=4,memory.G=30)
write.csv(rc,"RC.csv")

# Load betaNTI results
bNTI<-read.table("bNTI.group_viridescens_uninfected.csv",header=T,sep=",", check.names = F, row.names=1)
rc<-read.table("RC.csv",header=T,sep=",", check.names = F, row.names=1)
# Influences of different processes

#match names in bNTI and RC
rcc=rc[match(rownames(bNTI),rownames(rc)),match(colnames(bNTI),colnames(rc))]
bNTI.v<-as.vector(as.dist(bNTI))
rc.v<-as.vector(as.dist(rcc))
id.selectna<-(bNTI.v<=2&bNTI.v>=(-2))
num.pair<-length(bNTI.v)
select.h<-sum(bNTI.v>2)/num.pair
select.l<-sum(bNTI.v<(-2))/num.pair
disper.h<-sum(rc.v[id.selectna]>0.95)/num.pair
disper.l<-sum(rc.v[id.selectna]<(-0.95))/num.pair
drift<-sum(rc.v[id.selectna]<=0.95&rc.v[id.selectna]>=(-0.95))/num.pair
res=data.frame(select.h,select.l,disper.h,disper.l,drift,num.pair)
write.csv(res,"Processes.csv")

# select.h+select.l= Environment selection
select.h + select.l # 0.18
# disper.h= Dispersal limitation
disper.h #0.13
# disper.l= Homogenizing dispersal
disper.l #0.036
# drift= Ecological drift
drift #0.64
# num.pair= number of paired comparisons
