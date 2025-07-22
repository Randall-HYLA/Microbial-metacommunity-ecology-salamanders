# COMMUNITY ASSEMBLY ANALYSIS

# Part 1: Calculate bNTI

#Load files

#Load ASV table as matrix
otu.table <- read.csv("forest_otu_table_rare.csv", header = T, check.names = F, row.names = 1)
otu.matrix <- as.matrix(otu.table)

samp.meta <- read.csv("forest_metadata_rare.csv", head = T, sep = ",")

#Load phylogenetic tree
library(picante)
tree <- read.tree(file="forest_final_root.nwk")
tree <- ape::multi2di(tree)

# Calculate distance
dis <- cophenetic.phylo(tree)# calculate phylogenetic distance matrix
spname <- rownames(dis)# names of species with phylogenetic distance
colnames(dis) <- rownames(dis)

# Match the data
comm <- otu.matrix[match(spname,rownames(otu.matrix)), 1:ncol(otu.matrix)]# remove the species without phylogenetic information.
rownames(comm) <- spname
comm <- comm[,match(samp.meta[,1],colnames(comm))]# sort community data by sample grouping file. 
comm <- t(comm) # need transposition

# load function of betaNTI, for Hydra change this to full code
bNTI<-function(comm, dis, samp.group=NA, weighted=c(TRUE,FALSE),grouping=c(FALSE,TRUE),rand=1000,output.bMNTD=c(FALSE,TRUE))
{
  library(picante)
  samp.name=rownames(comm)
  result=data.frame(matrix(NA,length(samp.name),length(samp.name)))
  colnames(result)=samp.name
  rownames(result)=samp.name
  result.mntd=result
  
  if(grouping)
  {
    # calculate betaMNTD and betaNTI within group #
    group=levels(as.factor(samp.group[,2]))
    group.n=length(group)
    bMNTD=list()
    bNTI=list()
    pbar <- txtProgressBar(min = 0, max = 20, style = 3)# progress bar for large data
    for(m in 1:group.n)
    {
      samp.namex=samp.group[samp.group[,2]==group[m],1]#choose group
      if(length(samp.namex)<2){bNTI[[m]]=NA;bMNTD[[m]]=NA}else{
        comx=comm[match(samp.namex,rownames(comm)),]#remove others
        comx=comx[,colSums(comx != 0, na.rm = TRUE) > 0]#remove undetected species in this group, this part about remove NA, I added it myself
        spnamex=colnames(comx)
        disx=dis[match(spnamex,rownames(dis)),match(spnamex,colnames(dis))]#remove undetected species
        bMNTD.obs<-comdistnt(comx, disx, abundance.weighted = weighted, exclude.conspecifics = FALSE)# calculate observed betaMNTD.
        bMNTD[[m]]=as.matrix(bMNTD.obs)
        bMNTD.rand=array(dim=c(length(samp.namex),length(samp.namex),rand))
        for(i in 1:rand)
        {
          rand.namex=sample(spnamex)
          #disx.rand=dis[match(rand.namex,spnamex),match(rand.namex,spnamex)]
          disx.rand=disx
          colnames(disx.rand)=rand.namex
          rownames(disx.rand)=rand.namex
          bMNTD.rand[,,i]<-as.matrix(comdistnt(comx, disx.rand, abundance.weighted = weighted, exclude.conspecifics = FALSE))
          setTxtProgressBar(pbar, round(((m-1)*rand+i)*20/(group.n*rand),0))
        }
        bNTI[[m]]=(bMNTD[[m]]-apply(bMNTD.rand,c(1,2),mean))/(apply(bMNTD.rand,c(1,2),sd))
      }
      result[match(rownames(bNTI[[m]]),samp.name),match(colnames(bNTI[[m]]),samp.name)]=bNTI[[m]]
      result.mntd[match(samp.namex,samp.name),match(samp.namex,samp.name)]=bMNTD[[m]]
    }
    close(pbar)
  }else{
    # calculate across all samples #
    bMNTD.obs<-as.matrix(comdistnt(comm, dis, abundance.weighted = weighted, exclude.conspecifics = FALSE)) # calculate observed betaMNTD.
    spname=colnames(comm)
    bMNTD.rand=array(dim=c(length(samp.name),length(samp.name),rand))
    pbar <- txtProgressBar(min = 0, max = 20, style = 3)# progress bar for large data
    for(i in 1:rand)
    {
      rand.name=sample(spname)
      dis.rand=dis
      colnames(dis.rand)=rand.name
      rownames(dis.rand)=rand.name
      bMNTD.rand[,,i]<-as.matrix(comdistnt(comm, dis.rand, abundance.weighted = weighted, exclude.conspecifics = FALSE))
      setTxtProgressBar(pbar, round((i*20/rand),0))
    }
    close(pbar)
    bNTI=(bMNTD.obs-apply(bMNTD.rand,c(1,2),mean))/(apply(bMNTD.rand,c(1,2),sd))
    result=bNTI
    result.mntd=bMNTD.obs
  }
  if(output.bMNTD)
  {
    output=list(betaNTI=result,betaMNTD=result.mntd)
  }else{
    output=result}
  output
}

# Observed betaMNTD across all samples
bMNTD <- picante::comdistnt(comm, dis, abundance.weighted = TRUE, exclude.conspecifics = FALSE) # calculate observed betaMNTD. pay attention to weight by abundance or not.
write.csv(as.matrix(bMNTD),file="bMNTD.csv")# output observed betaMNTD

# Calculate betaNTI across all samples
bNTI.all <- bNTI(comm, dis, samp.group=NA, weighted=TRUE,grouping=FALSE,rand=1000,output.bMNTD=FALSE)
write.csv(bNTI.all,file="bNTIall.csv")# output

# Calculate betaMNTD and betaNTI within group
bNTI.group <- bNTI(comm, dis, samp.group=samp.meta, weighted=TRUE,grouping=TRUE,rand=1000,output.bMNTD=TRUE)

write.csv(bNTI.group$betaNTI,file="bNTI.group_forest.csv")
write.csv(bNTI.group$betaMNTD,file="bMNTD.group_forest.csv")
