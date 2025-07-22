library(stringr)
library(dplyr)
library(phyloseq)
library(plyr)
library(tidyverse)
library(tidyr)
library(igraph)
library(bipartite)

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
(viridescens_bd <- full_dataset %>%
  subset_samples(Species == "viridescens") %>%
    subset_samples(Bd == "yes") %>%  
  subset_samples(Site_Name != "hone quarry"))

(viridescens_nobd <- full_dataset %>%
    subset_samples(Species == "viridescens") %>%
    subset_samples(Bd == "no") %>%  
    subset_samples(Site_Name != "hone quarry"))

(pond <-  full_dataset %>%
  subset_samples(Species == "pond") %>%
  subset_samples(Site_Name != "hone quarry"))

#Filter by most abundant 100 taxa

#Viridescens Bd+
asv_table_viri_bd <- as.data.frame(otu_table(viridescens_bd))
asv_abundance_viri_bd <- rowSums(asv_table_viri_bd)
sorted_asv_abundance_viri_bd <- sort(asv_abundance_viri_bd, decreasing = TRUE)
top_100_asvs_viri_bd <- names(sorted_asv_abundance_viri_bd)[1:100]
(viridescens_bd <- prune_taxa(top_100_asvs_viri_bd, viridescens_bd))
rm(asv_table_viri_bd,asv_abundance_viri_bd,sorted_asv_abundance_viri_bd,top_100_asvs_viri_bd)

#Viridescens Bd-
asv_table_viri_nobd <- as.data.frame(otu_table(viridescens_nobd))
asv_abundance_viri_nobd <- rowSums(asv_table_viri_nobd)
sorted_asv_abundance_viri_nobd <- sort(asv_abundance_viri_nobd, decreasing = TRUE)
top_100_asvs_viri_nobd <- names(sorted_asv_abundance_viri_nobd)[1:100]
(viridescens_nobd <- prune_taxa(top_100_asvs_viri_nobd, viridescens_nobd))
rm(asv_table_viri_nobd,asv_abundance_viri_nobd,sorted_asv_abundance_viri_nobd,top_100_asvs_viri_nobd)

#Pond
asv_table_pond <- as.data.frame(otu_table(pond))
asv_abundance_pond <- rowSums(asv_table_pond)
sorted_asv_abundance_pond <- sort(asv_abundance_pond, decreasing = TRUE)
top_100_asvs_pond <- names(sorted_asv_abundance_pond)[1:100]
(pond <- prune_taxa(top_100_asvs_pond, pond))
rm(asv_table_pond, asv_abundance_pond, sorted_asv_abundance_pond, top_100_asvs_pond)

# Merge the phyloseq objects
(merged_phyloseq <- merge_phyloseq(viridescens_bd, viridescens_nobd, pond, method = "inner")) 

#Create new variable
dir.create("networks_bipartite")

otu.viridescens.pond <- as.data.frame((as.table(otu_table(merged_phyloseq),"matrix"))) #ASV_table

otu.viridescens.pond <- otu.viridescens.pond %>%
  dplyr::rename(ASV = Var1, Sample = Var2)

write.csv(otu.viridescens.pond, "viridescens_otu_table.csv", row.names = T)

#Extract taxonomy from merged 
tax.viridescens.pond <- as(tax_table(merged_phyloseq),"matrix") #taxonomy_metadata
tax.viridescens.pond <- as.data.frame(tax.viridescens.pond)
head(tax.viridescens.pond)

write.csv(tax.viridescens.pond, "viridescens_tax_table.csv", row.names = T)

#Add taxonomy to otu/asv
viri_bd <- left_join(otu.viridescens.pond, tax.viridescens.pond, by = "ASV")

#Merge OTU info / Taxonomy Info + Metadata

metadata <- read.csv("Salamanders_metadata_final_2020_DEC13.csv")
getwd()

#Combine two or more columns / string
metadata$treatment <- paste(metadata$Species, "-", metadata$Bd)

metadata_viridescens <- metadata %>%
  dplyr::select(Sample, treatment)

viri_bd_pond <- left_join(viri_bd, metadata_viridescens, by = "Sample")

#Create cov variable / the mean freq of ASV per treatment
dfs <- ddply(viri_bd_pond, .(ASV, treatment), summarise, cov = mean(Freq))
dfs <- t(spread(dfs, treatment, cov,  drop=TRUE , fill = 0))
in.mat = t(apply(dfs[2:4,], 1, as.numeric)) #Numbers change depending on number of variables
colnames(in.mat) <- dfs[1,] #Incidence matrix for generating bipartite networks

# Incidence matrix results -------------------------
in.mat #217

#------------------------ IGRAPH ------------------------
#---- Generating graph object  ----- 
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

#

#------------------------ BIPARTITE ------------------------ 
#library(bipartite)
# Plot bipartite network using bipartite package
# data in incidence matrix format
plotweb(sortweb(in.mat, sort.order="inc"), method="normal")

# Plot in matrix format

#I havent been able to make the visweb look nice, but this is not needed
visweb(sortweb(in.mat, sort.order="dec"), type= "none", # change to nested or diagonal 
       labsize= 3, square= "interaction", text= "none", textsize= 4)

visweb(in.mat, plotsize=5)

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
(l.bc <- node.metrics$`lower level`[2]) # species strength 

## Network metrics
network.metrics <-  networklevel(round(in.mat))
# network.metrics # Want to know about the metrics? Call ?networklevel

# Exploring by metric / This takes time to process
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

#https://chrischizinski.github.io/rstats/igraph-ggplotll/

#visNetwork creates a dynamic network.

# Create the interactive plot using visNetwork
visNetwork(nodes, edges) %>%
  visIgraphLayout(layout = 'layout.davidson.harel') %>%
  visEdges(width = "width")

#Plot with ggplot2

#The following will generate the final network, to achieve this I will extract values from previous objects.

fr.all <- layout_with_dh(g) # From line 120 
fr.all.df <- as.data.frame(fr.all)  ## convert the layout to a data.frame
fr.all.df$species <- V(bipartite_graph)$name  ## add in the bacterial species codes

#Add taxonomy
tax.viridescens.pond <- tax.viridescens.pond %>%
  select(-Species) %>%
  dplyr::rename(species = ASV)

fr.all.df <- left_join(fr.all.df, tax.viridescens.pond, 
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

edges <- edges %>%
  mutate(from = fct_recode(from,
                              "Pond" = "pond - NA",
                              "Newts Bd+" = "viridescens - no",
                              "Newts Bd-" = "viridescens - yes"
  ))

fr.all.df <- fr.all.df %>%
  mutate(species = fct_recode(species,
                              "Pond" = "pond - NA",
                              "Newts Bd+" = "viridescens - no",
                              "Newts Bd-" = "viridescens - yes"
  ))

#Generate subsets for adjusting asv and type var
type_var <- fr.all.df[1:3, ] #Three threatments N. viridescens Bd+, N. viridescens Bd- and pond 
type_asv <- fr.all.df[4:220, ] #220 ASV

vector_list <- edges %>%
    filter(from == "Pond") %>%
  pull(to)

edges_newt_env_list <- edges %>%
  filter(from != "Pond") %>%
  filter(to %in% vector_list) %>%
  pull(to)

edges_newt_env <- edges %>%
  filter(to %in% edges_newt_env_list)

# Define the custom color palette for the points
color_palette <- c("#d62728","#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "lightblue","#CC79A7", "#D55E00", "#FF7F0E")

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
                     color = width)
                     ) +
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
ggsave("NEWTS_NETWORK.pdf", g1, width = 300,height = 200, units="mm",dpi = 600)
