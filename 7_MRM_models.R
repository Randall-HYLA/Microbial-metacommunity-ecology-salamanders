# Author Owen G. Osborne
library(phyloseq)
library(raster) # for getData
library(geodist) # for geodist
library(ade4) # for dudi.pca
library(usedist)

# source functions
devtools::source_url("https://raw.githubusercontent.com/ogosborne/salamander_phylosymbiosis/main/code/phylosymbiosis_funcs.R")

#### microbiome data
load("data/phyloseq_full_dataset_01prev_rarefy.RData")

### sals data
sal <- subset_samples(full_datasetRare, Sample_Type == "Salamander")
## add metadata
df <- data.frame(sample_data(sal))
# group (spp x Bd)
df$group <- paste(df$Genus, df$Species, c(yes="BdPos",no="Bdneg")[df$Bd], sep = ".")
# area-habitat combo 
df$area_hab <- paste(df$Area, df$Habitat, sep = ".")
# add to ps obj
df <- sample_data(df)
sample_data(sal) <- df
rm(df)
# keep only groups of interest
salGoI <- c("Eurycea.bislineata.Bdneg","Nothopthalmus.viridescens.Bdneg", "Nothopthalmus.viridescens.BdPos", "Plethodon.cinereus.Bdneg")

sal <- subset_samples(sal, group %in% salGoI)
# remove forest N. viridescens
sal <- prune_samples(x = sal, samples = !(sample_names(sal) %in% paste0("HOQ",15:19)))
# remove empty taxa
sal <- filter_taxa(sal, function(x) sum(x) != 0, TRUE)
# env only
env <- subset_samples(full_datasetRare, Sample_Type == "Environment")
df <- data.frame(sample_data(env))
# group (habitat)
df$group <- df$Habitat
# area-habitat combo 
df$area_hab <- paste(df$Area, df$Habitat, sep = ".")
# add to ps obj
df <- sample_data(df)
sample_data(env) <- df
rm(df)
# remove empty taxa
env <- filter_taxa(env, function(x) sum(x) != 0, TRUE)
#### whole ds microbiome distances
dists <- c("jaccard","bray")
dismats <- list()
DS <- c("env","sal")
for(D in DS){
  dismats[[paste(D, "bray", sep= ".")]] <- phyloseq::distance(get(D), method = "bray")
  dismats[[paste(D, "jaccard", sep= ".")]] <- phyloseq::distance(get(D), method = "jaccard", binary = TRUE)
}
### geographic distance
## metadata for filtered datasets 
sal_metadata <- data.frame(sample_data(sal))
env_metadata <- data.frame(sample_data(env))
all_metadata <- rbind(sal_metadata, env_metadata)
coords <- setNames(all_metadata[,c("longitude", "latitude")], c("x", "y"))
## get geodist
geog_dist <- geodist(coords, measure = "vincenty")
dimnames(geog_dist) <- list(rownames(all_metadata),rownames(all_metadata))
### climatic distance
# get bioclim data
bioclim <- getData("worldclim", var="bio", res = 0.5, lon = coords$x, lat = coords$y)
# extract values for sampling sites
bc_vals <- extract(bioclim, coords)
rownames(bc_vals) <- rownames(coords)
bc_vals <- as.data.frame(bc_vals)
rm(bioclim)
# add elevation
bc_vals$elev <- all_metadata$elevation
# run pca
bc_pca <- dudi.pca(bc_vals, scannf = FALSE, nf = 3)
# get climate distance matrix
clim_dist <- dist(bc_pca$li)
### environmental microbiome distance for each salamander sample
envm_dist <- list()
dists <- c("jaccard","bray")

for(d in dists){
  # first get the mean distance between environmental samples in each locality-habitat combination
  ds <- paste("env", d, sep = ".")
  my.envm <- summarise_mat(x = as.matrix(dismats[[ds]]), metadata = env_metadata, col = "area_hab", fun = "mean")
  # then set the diagonal to 0 so there is 0 distance between samples from the same locality and habitat. While samples vary within a habitat-locality, there is no 1-to-1 correspondence between environmental and salamander samples so there would otherwise just be the same average between all members of the same group.
  diag(my.envm) <- 0
  # then expand the matrix to give mean environmental microbiome distance between each pair of salamander samples (dismats$sal.jaccard is just used to get salamander sample names in order, distances are taken from my.envm)
  my.envm <- expand_mat(as.matrix(dismats$sal.jaccard), my.envm, metadata = sal_metadata, col = "area_hab")
  # convert to dist object
  my.envm <- as.dist(my.envm)
  # add to output list
  envm_dist[[paste(d, sep = ".")]] <- my.envm
}
# get data for each subset
salGoI <- c("Eurycea.bislineata.Bdneg","Nothopthalmus.viridescens.Bdneg", "Nothopthalmus.viridescens.BdPos", "Plethodon.cinereus.Bdneg")
envGoI <- c("forest", "pond", "stream")
dists <- c("jaccard","bray")
MRM_dat <- list()
# env subsets
for(i in envGoI){
  my.list <- list()
  # microbiome distance
  set <- row.names(env_metadata)[which(env_metadata$group == i)]
  for(d in dists){
    ds <- paste("env", d, sep = ".")
    my.list[[d]] <- dist_subset(dismats[[ds]], set)
  }
  # climate and geog distance
  my.list$clim_dist <-dist_subset(clim_dist, set) 
  my.list$geog_dist <-dist_subset(geog_dist, set) 
  # scale all dismats
  my.list <-  lapply( my.list, function(x) scales::rescale(as.dist(x), to = 0:1))
  # add to output list
  MRM_dat[[paste0("env.",i)]] <- my.list
}
# sal subsets
for(i in salGoI){
  my.list <- list()
  # microbiome distance
  set <- row.names(sal_metadata)[which(sal_metadata$group == i)]
  for(d in dists){
    ds <- paste("sal", d, sep = ".")
    my.list[[d]] <- dist_subset(dismats[[ds]], set)
  }
  # environmental microbiome distance
  for(d in dists){
    ds.in <- paste(d, sep = ".")
    ds.out <- paste("envm", d, sep = ".")
    my.list[[ds.out]] <- dist_subset(envm_dist[[ds.in]], set) 
  }
  # climate and geog distance
  my.list$clim_dist <-dist_subset(clim_dist, set) 
  my.list$geog_dist <-dist_subset(geog_dist, set) 
  # standardise all dismats
  my.list <- lapply(my.list, function(x) scales::rescale(as.dist(x), to = 0:1))#function(x) (x - mean(x))/sd(x) )
  # add to output list
  MRM_dat[[paste0("sal.",i)]] <- my.list
}
# check all covariates are in the same order
checks <- c()
for(setN in names(MRM_dat)){
  varNs <- names(MRM_dat[[setN]])
  for (varN in varNs){
   checks <- c(checks, all(rownames(as.matrix(MRM_dat[[setN]][[1]])) == rownames(as.matrix(MRM_dat[[setN]][[varN]]))))
   checks <- c(checks, all(colnames(as.matrix(MRM_dat[[setN]][[1]])) == colnames(as.matrix(MRM_dat[[setN]][[varN]]))))
  }
}
all(checks)
# save
dir.create("results")
save(dismats, file = "results/dismats2.RData")
save(MRM_dat, file = "results/MRM_dat2.RData")
save(list = c("sal", "env", "sal_metadata", "env_metadata"), file = "results/subset_ps2.RData")
### run MRM
# clear environment
invisible(lapply(paste0('package:', names(sessionInfo()$otherPkgs)), detach, character.only=TRUE, unload=TRUE, force=TRUE))
rm(list=ls())
library(ecodist)
library(gtools)
library(RColorBrewer)
load("results/MRM_dat2.RData")
# run 
salGoI <- c("Eurycea.bislineata.Bdneg","Nothopthalmus.viridescens.Bdneg", "Nothopthalmus.viridescens.BdPos", "Plethodon.cinereus.Bdneg")
envGoI <- c("forest", "pond", "stream")
dists <- c("jaccard","bray")
set.seed(878703)
# env
MRM_res_env <- list()
for(i in envGoI){
  for(d in dists){
    ds <- paste0("env.",i)
    my.dat <- MRM_dat[[ds]]
    my.dist <- my.dat[[d]]
    clim <- my.dat$clim_dist
    geog <- my.dat$geog_dist
    my.res <- MRM(my.dist ~ geog + clim, nperm = 10000, mrank = TRUE)
    my.df <- data.frame(group = i,
                        distance = d,
                        R2_perm = my.res$r.squared["R2"],
                        Pval_perm = my.res$r.squared["pval"],
                        F.stat = my.res$F.test[["F"]],
                        F.pval = my.res$F.test[["F.pval"]],
                        clim.coeff = my.res$coef["clim", 1],
                        clim.pval = my.res$coef["clim", 2],
                        geog.coeff = my.res$coef["geog", 1],
                        geog.pval = my.res$coef["geog", 2]
                        )
    MRM_res_env[[paste(i,d,sep = ".")]] <-  my.df
  }
}
MRM_res_env <- do.call(rbind, MRM_res_env)
# sal
MRM_res_sal <- list()
for(i in salGoI){
  for(d in dists){
    for(D in dists){
      ds <- paste0("sal.",i)
      my.dat <- MRM_dat[[ds]]
      my.dist <- my.dat[[d]]
      clim <- my.dat$clim_dist
      geog <- my.dat$geog_dist
      envm <- my.dat[[paste("envm", D, sep = ".")]]
      my.res <- MRM(my.dist ~ geog + clim + envm, nperm = 10000, mrank = TRUE)
      my.df <- data.frame(group = i,
                          saldistance = d,
                          envdistance = D,
                          R2_perm = my.res$r.squared["R2"],
                          Pval_perm = my.res$r.squared["pval"],
                          F.stat = my.res$F.test[["F"]],
                          F.pval = my.res$F.test[["F.pval"]],
                          clim.coeff = my.res$coef["clim", 1],
                          clim.pval = my.res$coef["clim", 2],
                          geog.coeff = my.res$coef["geog", 1],
                          geog.pval = my.res$coef["geog", 2],
                          envm.coeff = my.res$coef["envm", 1],
                          envm.pval = my.res$coef["envm", 2]
      )
      MRM_res_sal[[paste0(i, "_saldist.", d, "_envdist.", D)]] <-  my.df 
    }
  }
}
MRM_res_sal <- do.call(rbind, MRM_res_sal)
# save results
write.csv(MRM_res_env, file = "results/MRM_res_env.csv", row.names = F, quote = F)
write.csv(MRM_res_sal, file = "results/MRM_res_sal_all_metric_combos.csv", row.names = F, quote = F)
# subset
plotdist <- "bray"
#plotdist <- "jaccard"
MRM_res_env_ss <- MRM_res_env[which(MRM_res_env$distance == plotdist),]
MRM_res_sal_ss <- MRM_res_sal[which(MRM_res_sal$distance == plotdist),]
# get ranges
range.env <- scales::expand_range(range(c(MRM_res_env_ss$clim.coeff, MRM_res_env_ss$geog.coeff)), mul =0.1)
range.sal <- scales::expand_range(range(c(MRM_res_sal_ss$clim.coeff, MRM_res_sal_ss$geog.coeff, MRM_res_sal_ss$envm.coeff)), mul =0.1)
# stars function
dostars <- function(bp, dat, ns = F){
  ypos  <- c(dat[,1:3])
  ypos <- ifelse(ypos > 0, ypos + 0.03, 0.03)
  stars <- stars.pval(dat[,4:6])
  stars[which(stars == ".")] <- " "
  if(ns) stars[which(stars == " ")] <- "n.s"
  text(x = bp, y = ypos, labels = stars, cex = 2)
}
#### Plot coefficients
# colour
mypal <- c("#4CAE4A", # green = environmental selection = clim
           "#E4191A", # red = dispersal limitation = geog
           "#2E78B5")# blue = ecological drift = envm
# layout
layout(matrix(1:6,byrow=T,ncol=3))
par(mar=c(3, 6.1, 2, 1))
# forest
# data
dat <- MRM_res_env_ss[which(MRM_res_env_ss$group == "forest"), c("clim.coeff", "geog.coeff","clim.pval","geog.pval")]
dat <- as.matrix(data.frame(clim.coeff = dat[,1], geog.coeff = dat[,2], envm.coeff = 0, clim.pval = dat[,3], geog.pval = dat[,4], envm.pval = 1))
# barplot
bp <- barplot(dat[,1:3],  beside = T, names.arg = c("", "", ""), ylab = "Regression coefficient", yaxt = "s", col = mypal, ylim =  range.env, main = "Forest", cex.lab = 2, cex.main = 2, cex.names = 2, font.main = 1)
abline(h=0, lwd = 1.5)
dostars(bp, dat)
# stream
# data
dat <- MRM_res_env_ss[which(MRM_res_env_ss$group == "stream"), c("clim.coeff", "geog.coeff","clim.pval","geog.pval")]
dat <- as.matrix(data.frame(clim.coeff = dat[,1], geog.coeff = dat[,2], envm.coeff = 0, clim.pval = dat[,3], geog.pval = dat[,4], envm.pval = 1))
# barplot
bp <- barplot(dat[,1:3],  beside = T, names.arg = c("", "", ""), ylab = "", yaxt = "n", col = mypal, ylim = range.env, main = "Stream", cex.main = 2, cex.names = 2, font.main = 1)
abline(h=0, lwd = 1.5)
dostars(bp, dat)
# pond
# data
dat <- MRM_res_env_ss[which(MRM_res_env_ss$group == "pond"), c("clim.coeff", "geog.coeff","clim.pval","geog.pval")]
dat <- as.matrix(data.frame(clim.coeff = dat[,1], geog.coeff = dat[,2], envm.coeff = 0, clim.pval = dat[,3], geog.pval = dat[,4], envm.pval = 1))
# barplot
bp <- barplot(dat[,1:3],  beside = T, names.arg = c("", "", ""), ylab = "", yaxt = "n", col = mypal, ylim = range.env, main = "Pond", cex.main = 2, cex.names = 2, font.main = 1)
abline(h=0, lwd = 1.5)
dostars(bp, dat)
# P. cin
# data
dat <- as.matrix(MRM_res_sal_ss[which(MRM_res_sal_ss$group == "Plethodon.cinereus.Bdneg"), c("clim.coeff", "geog.coeff", "envm.coeff", "clim.pval", "geog.pval", "envm.pval")])
# barplot
bp <- barplot(dat[,1:3],  beside = T, names.arg = c("Clim", "Geog", "Envm"), las = 1, ylab = "Regression coefficient", yaxt = "s", col = mypal, cex.lab = 2, ylim = range.sal, main = expression(paste(italic("P. cinereus"), "  Bd-")), cex.main = 2, cex.names = 2)
abline(h=0, lwd = 1.5)
dostars(bp, dat)
# E. bis
# data
dat <- as.matrix(MRM_res_sal_ss[which(MRM_res_sal_ss$group == "Eurycea.bislineata.Bdneg"), c("clim.coeff", "geog.coeff", "envm.coeff", "clim.pval", "geog.pval", "envm.pval")])
# barplot
bp <- barplot(dat[,1:3],  beside = T, names.arg =  c("Clim", "Geog", "Envm"), las = 1, ylab = "", yaxt = "n", col = mypal, cex.lab = 1.5, ylim = range.sal, main = expression(paste(italic("E. bislineata"), "  Bd-")), cex.main = 2, cex.names = 2)
abline(h=0, lwd = 1.5)
dostars(bp, dat)
# N. vir
# data
dat1 <- as.matrix(MRM_res_sal_ss[which(MRM_res_sal_ss$group == "Nothopthalmus.viridescens.BdPos" ), c("clim.coeff", "geog.coeff", "envm.coeff", "clim.pval", "geog.pval", "envm.pval")])
dat2 <- as.matrix(MRM_res_sal_ss[which(MRM_res_sal_ss$group == "Nothopthalmus.viridescens.Bdneg" ), c("clim.coeff", "geog.coeff", "envm.coeff", "clim.pval", "geog.pval", "envm.pval")])
dat <- rbind(dat1,dat2)
# barplot
bp <- barplot(dat[,1:3],  beside = T, names.arg = c("Clim", "Geog", "Envm"), las = 1, ylab = "", yaxt = "n", col = mypal[c(1,1,2,2,3,3)], density = c(20,1000), angle = 45, cex.lab = 1.5, ylim = range.sal, main = expression(italic("N. viridescens")), cex.main = 2, cex.names = 2)
abline(h=0, lwd = 1.5)
dostars(bp, dat)
legend("topleft", legend = c("Bd+", "Bd-"), fill = "black", bty = "n", cex = 2.2, density = c(20,NA), angle = 45)
# Jaccard Bray correlation plots
par(mfrow = c(4,2), mar = c(5.1,4.1,4.1,2.1))
plot(MRM_dat$sal.Eurycea.bislineata.Bdneg$jaccard, MRM_dat$sal.Eurycea.bislineata.Bdneg$bray, xlab = "Jaccard", ylab = "Bray", main = "E.bis")
plot(MRM_dat$sal.Plethodon.cinereus.Bdneg$jaccard, MRM_dat$sal.Plethodon.cinereus.Bdneg$bray, xlab = "Jaccard", ylab = "Bray", main = "P.cin")
plot(MRM_dat$sal.Nothopthalmus.viridescens.Bdneg$jaccard, MRM_dat$sal.Nothopthalmus.viridescens.Bdneg$bray, xlab = "Jaccard", ylab = "Bray", main = "N.vir Bd-")
plot(MRM_dat$sal.Nothopthalmus.viridescens.BdPos$jaccard, MRM_dat$sal.Nothopthalmus.viridescens.BdPos$bray, xlab = "Jaccard", ylab = "Bray", main = "N.vir Bd+")
plot(MRM_dat$env.forest$jaccard, MRM_dat$env.forest$bray, xlab = "Jaccard", ylab = "Bray", main = "Forest")
plot(MRM_dat$env.pond$jaccard, MRM_dat$env.pond$bray, xlab = "Jaccard", ylab = "Bray", main = "Pond")
plot(MRM_dat$env.stream$jaccard, MRM_dat$env.stream$bray, xlab = "Jaccard", ylab = "Bray", main = "Stream")
