library(vegan);library(parallel);library(microbiome);library(picante);library(ape)

# Part 2: 
com <- read.csv("pond_otu_table_rare.csv", header = T, check.names = F, row.names = 1)
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
bNTI<-read.table("bNTI.group_pond.csv",header=T,sep=",", check.names = F, row.names=1)
rc <- read.table("RC.csv",header=T,sep=",", check.names = F, row.names=1)
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
select.h + select.l # 0.52
# disper.h= Dispersal limitation
disper.h #0.34
# disper.l= Homogenizing dispersal
disper.l #0.01
# drift= Ecological drift
drift #0.13
# num.pair= number of paired comparisons
