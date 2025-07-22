library(stringr)
library(dplyr)
library(phyloseq)
library(plyr)
library(tidyverse)
library(tidyr)

load("phyloseq_full_dataset_01prev_rarefy.RData") #Unrarefied dataset

full_dataset <- full_datasetRare
rm(full_datasetRare)

sort(sample_sums(full_dataset))
max(sample_sums(full_dataset))/min(sample_sums(full_dataset))

## Rename UNCLASSIFIED taxa

#Rename "NA" to "unspecified + [last identified taxa]"
tax.clean <- data.frame(tax_table(full_dataset)) 
for (i in 1:7){ tax.clean[,i] <- as.character(tax.clean[,i])} 
tax.clean[is.na(tax.clean)] <- "Unassigned" #Change NA to and empty string
tax_table(full_dataset) <- as.matrix(tax.clean) 

#Create subsets for filtering
(bislineata <- full_dataset %>%
  subset_samples(Species == "bislineata" & Bd == "no"))
  
(stream <-  full_dataset %>%
  subset_samples(Species == "stream"))

#Filter by most abundant 100 taxa

#bislineata
asv_table_bislineata <- as.data.frame(otu_table(bislineata))
asv_abundance_bislineata <- rowSums(asv_table_bislineata)
sorted_asv_abundance_bislineata <- sort(asv_abundance_bislineata, decreasing = TRUE)
top_100_asvs_bislineata <- names(sorted_asv_abundance_bislineata)[1:100]
(bislineata <- prune_taxa(top_100_asvs_bislineata, bislineata))
rm(asv_table_bislineata,asv_abundance_bislineata,sorted_asv_abundance_bislineata,top_100_asvs_bislineata)

#stream
asv_table_stream <- as.data.frame(otu_table(stream))
asv_abundance_stream <- rowSums(asv_table_stream)
sorted_asv_abundance_stream <- sort(asv_abundance_stream, decreasing = TRUE)
top_100_asvs_stream <- names(sorted_asv_abundance_stream)[1:100]
(stream <- prune_taxa(top_100_asvs_stream, stream))
rm(asv_table_stream, asv_abundance_stream, sorted_asv_abundance_stream, top_100_asvs_stream)

# Merge the phyloseq objects
(merged_phyloseq <- merge_phyloseq(bislineata, stream, method = "inner")) 

#Create new variable
dir.create("networks_bipartite")

otu.bislineata.stream <- as.data.frame((as.table(otu_table(merged_phyloseq),"matrix"))) #ASV_table

otu.bislineata.stream <- otu.bislineata.stream %>%
  dplyr::rename(ASV = Var1, Sample = Var2)

write.csv(otu.bislineata.stream, "bislineata_otu_table.csv", row.names = T)

#Extract taxonomy from merged 
tax.bislineata.stream <- as(tax_table(merged_phyloseq),"matrix") #taxonomy_metadata
tax.bislineata.stream <- as.data.frame(tax.bislineata.stream)
head(tax.bislineata.stream)

write.csv(tax.bislineata.stream, "bislineata_tax_table.csv", row.names = T)

bislineata_bd <- left_join(otu.bislineata.stream, tax.bislineata.stream, by = "ASV")

#Merge OTU info / Taxonomy Info + Metadata
metadata <- read.csv("Salamanders_metadata_final_2020_DEC13.csv")

#Combine two or more columns / string
metadata$treatment <- paste(metadata$Species)

metadata_bislineata <- metadata %>%
  dplyr::select(Sample, treatment)

bislineata_bd_stream <- left_join(bislineata_bd, metadata_bislineata, by = "Sample")

#Create new var
dfs <- ddply(bislineata_bd_stream, .(ASV, treatment), summarise, cov = mean(Freq))
dfs <- t(spread(dfs, treatment, cov,  drop=TRUE , fill = 0))

in.mat = t(apply(dfs[2:3,], 1, as.numeric)) #Numbers change depending on number of variables
colnames(in.mat) <- dfs[1,]

# Incidence matrix results -------------------------
in.mat #181

#------------------------ IGRAPH ------------------------
#---- Generating graph object  ----- 
library(igraph)

g <- graph.incidence(in.mat, weighted= TRUE)

# Network attributes
V(g)$name # Check the vertex names 
V(g)$type # Check vertex types 

# plotting
plot(g, vertex.label.color='black', vertex.label.dist= ((V(g)$type * 4)-2)*-1  , layout = layout_as_bipartite)
plot(g, vertex.label.color='black')

#---- Adding attributes
# Shapes to nodes
shapes = c(rep("square", length(V(g)$type[V(g)$type == "FALSE"])), 
           rep("circle", length(V(g)$type[V(g)$type == "TRUE"])))
# plotting
plot(g, vertex.shape=shapes, vertex.label.color='black')

#------------------------
# Add colors to nodes

V(g)$color <- c(rep("#046A38", length(V(g)$type[V(g)$type == "FALSE"])), 
                rep("#FA4616", length(V(g)$type[V(g)$type == "TRUE"])))

# plotting
plot(g, vertex.shape=shapes, vertex.label.color='black')

#------------------------ BIPARTITE ------------------------ 
library(bipartite)
# Plot bipartite network using bipartite package
# data in incidence matrix format
plotweb(sortweb(in.mat, sort.order="inc"), method="normal")

# Plot in matrix format
visweb(sortweb(in.mat, sort.order="dec"), type= "none", # change to nested or diagonal 
       labsize= 3, square= "interaction", text= "none", textsize= 4)

visweb(in.mat)

#----- Calculating metrics ----

## Node metrics
node.metrics <- specieslevel(round(in.mat))

# Exploring metrics
str(node.metrics)
# How many levels are in the list?

# node.metrics$`higher level` # Want to know about the metrics? Call ?specieslevel

# Exploring $`higher level`
(h.nd <- node.metrics$`higher level`[1]) # node degree OR node.metrics$`higher level`$degree
(h.bc <- node.metrics$`higher level`[2]) # species strength 

# Exploring $`lower level`
(l.nd <- node.metrics$`lower level`[1]) # node degree
(l.bc <- node.metrics$`lower level`[2]) # habitat strength 

## Network metrics
network.metrics <-  networklevel(round(in.mat))
# network.metrics # Want to know about the metrics? Call ?networklevel

# Exploring by metric
network.metrics["connectance"] # Connectance
network.metrics["weighted nestedness"] # Nestedness *weighted

# Computing modularity
computeModules(in.mat) # Default method: Becket
(modularity <-  LPA_wb_plus(in.mat)) 

mod <- convert2moduleWeb(in.mat, modularity)

#----Plotting with attributes ----

# Combining node attributes from bipartite and igraph 
V(g)$xx <- c(unlist(l.nd), unlist(h.bc)+5) # adding node degree + species strength plus a constant

# Types of layout algorithms
#default: Kamada-Kawai 
plot(g, vertex.shape=shapes, edge.width= log(E(g)$weight), vertex.size= as.numeric(V(g)$xx))
# Davidson-Harel

plot(g, vertex.shape=shapes, edge.width= log(E(g)$weight), vertex.size= as.numeric(V(g)$xx), layout=layout_with_dh(g) )

 # Distributed Recursive
plot(g, vertex.shape=shapes, edge.width= log(E(g)$weight), vertex.size= as.numeric(V(g)$xx), layout=layout_with_drl(g) )

# Create the bipartite graph from the sparse matrix
bipartite_graph <- graph_from_incidence_matrix(in.mat, directed = FALSE)

library(visNetwork)
# Prepare data for visNetwork
nodes <- data.frame(id = V(bipartite_graph)$name, group = V(bipartite_graph)$type)
nodes$color.background <- ifelse(nodes$group, "lightblue", "red")
nodes$label <- nodes$id

edges <- get.data.frame(bipartite_graph, what = "edges")
edges$width <- log(E(g)$weight)
E(g)$weight
#https://chrischizinski.github.io/rstats/igraph-ggplotll/

# Create the interactive plot using visNetwork
visNetwork(nodes, edges) %>%
  visIgraphLayout(layout = 'layout.davidson.harel') %>%
  visEdges(width = "width")

#Plot with ggplot2
fr.all <- layout_with_dh(g)
fr.all.df <- as.data.frame(fr.all)  ## convert the layout to a data.frame
fr.all.df$species <- V(bipartite_graph)$name  ## add in the species codes

#Add taxonomy
tax.bislineata.stream <- tax.bislineata.stream %>%
  select(-Species) %>%
  dplyr::rename(species = ASV)

fr.all.df <- left_join(fr.all.df, tax.bislineata.stream, 
                       by = "species")

#Add type / FALSE = variable type - TRUE = bacteria
fr.all.df$type <- V(bipartite_graph)$type
fr.all.df  ## display the x (V1) and y (V2) coordinates for each of the nodes.

edges <- get.data.frame(bipartite_graph, what = "edges")
edges$width <- log(E(g)$weight) # get the edge information using the get.data.frame function

edges$from.x <- fr.all.df$V1[match(edges$from, fr.all.df$species)]  #  match the from locations from the node data.frame we previously connected
edges$from.y <- fr.all.df$V2[match(edges$from, fr.all.df$species)]
edges$to.x <- fr.all.df$V1[match(edges$to, fr.all.df$species)]  #  match the to locations from the node data.frame we previously connected
edges$to.y <- fr.all.df$V2[match(edges$to, fr.all.df$species)]

#Generate subsets for adjusting asv and type var
type_var <- fr.all.df[1:2, ]
type_asv <- fr.all.df[3:183, ]

vector_list <- edges %>%
    filter(from == "stream") %>%
  pull(to)

edges_newt_env_list <- edges %>%
  filter(from != "stream") %>%
  filter(to %in% vector_list) %>%
  pull(to)

edges_newt_env <- edges %>%
  filter(to %in% edges_newt_env_list)

color_palette <- c("#d62728", "#E69F00", "#56B4E9", "#0072B2", "#9467BD","#FF7F0E")


library(ggnewscale)


(g1 <- ggplot() +
  geom_segment(data=edges,
               aes(x=from.x, 
                   xend = to.x, 
                   y=from.y, 
                   yend = to.y, 
                   size = width),
                   color = "lightgray") +
    geom_segment(data=edges_newt_env,
                 aes(x=from.x, 
                     xend = to.x, 
                     y=from.y, 
                     yend = to.y, 
                     size = width,
                     color = width)) +
    scale_color_gradient(low = "yellow", high = "red") +
    scale_size_continuous(range = c(0.01, 1.0)) +
  
    new_scale_colour() + #REALLY IMPORTANT TO SET UP NEW COLOR SCALE
    
    geom_point(data=type_asv, 
             aes(x=V1,y=V2, color=Phylum), 
             size=3) +
    scale_color_manual(values = color_palette) +  # Set color palette
    geom_point(data=type_var, 
             aes(x=V1,y=V2),
             size=14,
             shape= 15,
             colour="azure4") +  # adds a black border around the nodes
  geom_point(data=type_var, 
             aes(x=V1,y=V2),
             color="gray",
             shape = 15,
             size=12) +
  #scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 15)) +
  geom_text(data=type_var, 
            aes(x=V1,y=V2,label=species), 
            size = 3) + # add the node labels
  scale_x_continuous(expand=c(0,1))+  # expand the x limits 
  scale_y_continuous(expand=c(0,1))+ # expand the y limits
  theme_bw() +  # use the ggplot black and white theme
  theme(
    axis.text.x = element_blank(),  # remove x-axis text
    axis.text.y = element_blank(), # remove y-axis text
    axis.ticks = element_blank(),  # remove axis ticks
    axis.title.x = element_blank(), # remove x-axis labels
    axis.title.y = element_blank(), # remove y-axis labels
    panel.background = element_blank(), 
    panel.border =element_blank(), 
    panel.grid.major = element_blank(),  #remove major-grid labels
    panel.grid.minor = element_blank(),  #remove minor-grid labels
    plot.background = element_blank())
)
  
# Export the plot using ggsave()
ggsave("BISLINEATA_NETWORK.pdf", g1, width = 300,height = 200, units="mm",dpi = 600)

#Create heatmap of ASV shared between salamanders and environment

#Get a list of ASVs shared 
shared_asvs <- unique(edges_newt_env$to)

data_heatmap <- ddply(bislineata_bd_stream, .(ASV, treatment), summarise, cov = mean(Freq))
data_heatmap <- subset(data_heatmap, ASV %in% shared_asvs)

data_heatmap <- data_heatmap %>% 
  filter(cov != 0)

data_heatmap$logfreq <- log10(data_heatmap$cov+1) 
data_heatmap$cov <- round(data_heatmap$cov, 1)

data_heatmap <- left_join(data_heatmap,tax.clean, by = "ASV")
data_heatmap$ASV_tax <- paste(data_heatmap$Class, data_heatmap$Order, data_heatmap$Family, data_heatmap$ASV)

# Reorder ASV_tax in descending alphabetical order
data_heatmap <- data_heatmap %>%
  mutate(ASV_tax = factor(ASV_tax, levels = rev(sort(unique(ASV_tax)))))

data_heatmap$treatment <- factor(data_heatmap$treatment, levels = c("stream", "bislineata"))

(plot_indicator <- ggplot(data_heatmap, aes(x = treatment, y = ASV_tax)) +
    geom_tile(aes(fill = logfreq), colour="white", size=0.30) +
    xlab("") +
    geom_text(aes(label=cov, size=2)) +
    #scale_fill_gradientn(colours = hm.palette(100))+
    scale_fill_gradient(low = "yellow1", high = "red", name="Log\nabundance")+
    theme_bw() +
    theme(axis.title.x = element_text(face="bold", color = "black", size=12),
          axis.title.y=element_blank(),
          axis.text.x=element_text(colour="black", size=10),
          legend.position = "right",
          axis.text.y=element_text(colour="black", size = 10),
          axis.line = element_line(colour = "black"),
          strip.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank())
)

ggsave("heatmap_bislineata_stream_taxa.pdf", plot_indicator, width = 250, height = 140, units="mm",dpi = 300)

# Fit a generalized mixed model
model <- lmer(width ~ from + (1 | to), data = 
                edges_newt_env)

# Summarize the model
summary(model)

library(car)
Anova(model, test.statistic = "F")
