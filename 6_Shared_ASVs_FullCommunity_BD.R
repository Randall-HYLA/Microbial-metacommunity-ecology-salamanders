# ASV overlap analysis 

library(stringr)
library(dplyr)
library(phyloseq)
library(ggplot2)

load("phyloseq_full_dataset_01prev_rarefy.RData")

# Filter data
full_datasetRare <- full_datasetRare %>%
  subset_samples(Species != "bislineata" | Bd != "yes")

full_datasetRare_filt  <- full_datasetRare  %>%
  subset_samples(Species %in% c("viridescens", "cinereus", "bislineata", "forest", "pond", "stream"))

full_datasetRare_filt@sam_data$species_new <- paste(full_datasetRare_filt@sam_data$Species, full_datasetRare_filt@sam_data$Bd)

full_datasetRare_filt <- prune_taxa(taxa_sums(full_datasetRare_filt) > 0, full_datasetRare_filt)

# Melt data
otu_glom <- tax_glom(full_datasetRare_filt, taxrank = "ASV")
otu_melt <- psmelt(otu_glom)

# Get OTU abundances by group
group_abundances <- function(group_label, species_condition) {
  otu_abund <- subset(otu_melt, {{species_condition}}) %>%
    select(Sample, Abundance, OTU) %>%
    group_by(OTU) %>%
    summarize(OTU_sum = sum(Abundance)) %>%
    filter(OTU_sum > 0)
  colnames(otu_abund)[2] <- paste0("OTU_sum_", group_label)
  otu_abund
}

otu_abund_viridescens_bd    <- group_abundances("viridescens_bd", species_new == "viridescens yes")
otu_abund_viridescens_nobd  <- group_abundances("viridescens_nobd", species_new == "viridescens no")
otu_abund_bislineata        <- group_abundances("bislineata", species_new == "bislineata no")
otu_abund_cinereus          <- group_abundances("cinereus", species_new == "cinereus no")
otu_abund_forest            <- group_abundances("forest", sample_Species == "forest")
otu_abund_pond              <- group_abundances("pond", sample_Species == "pond")
otu_abund_stream            <- group_abundances("stream", sample_Species == "stream")

# Pairwise OTU comparisons
otu_compare <- function(group1, group2, label1, label2, comparison_label) {
  only_1 <- anti_join(group1, group2, by = "OTU")
  only_2 <- anti_join(group2, group1, by = "OTU")
  both   <- inner_join(group1, group2, by = "OTU") %>%
    mutate(Comparison = comparison_label,
           Overlap = paste("Both", label1, "and", label2),
           TotalSum = across(starts_with("OTU_sum"), sum, .names = "TotalSum"))
  list(only_1 = only_1, only_2 = only_2, both = both)
}

comp1 <- otu_compare(otu_abund_viridescens_bd, otu_abund_pond, "NewtBd", "Pond", "N. viridescens Bd \n vs. \n Pond")
comp2 <- otu_compare(otu_abund_viridescens_nobd, otu_abund_pond, "NewtnoBd", "Pond", "N. viridescens noBd \n vs. \n Pond")
comp3 <- otu_compare(otu_abund_cinereus, otu_abund_forest, "Cinereus", "Forest", "P. cinereus \n vs. \n Forest")
comp4 <- otu_compare(otu_abund_bislineata, otu_abund_stream, "Bislineata", "Stream", "E. bislineata \n vs. \n Stream")

# Create OTU overlap table
otu_overlap_env <- data.frame(
  SampleType = c("NewtBd", "NewtBd", "Pond", "Pond", "NewtnoBd", "NewtnoBd", "Pond", "Pond", 
                 "Bislineata", "Bislineata", "Stream", "Stream", "cinereus", "cinereus", "Forest", "Forest"),
  Type = c("NewtBd only", "Shared", "Pond only", "Shared",
           "NewtnoBd only", "Shared", "Pond only", "Shared",
           "Bislineata only", "Shared", "Stream only", "Shared",
           "Cinereus only", "Shared", "Forest only", "Shared"),
  NumOTUs = c(nrow(comp1$only_1), nrow(comp1$both),
              nrow(comp1$only_2), nrow(comp1$both),
              nrow(comp2$only_1), nrow(comp2$both),
              nrow(comp2$only_2), nrow(comp2$both),
              nrow(comp4$only_1), nrow(comp4$both),
              nrow(comp4$only_2), nrow(comp4$both),
              nrow(comp3$only_1), nrow(comp3$both),
              nrow(comp3$only_2), nrow(comp3$both)),
  Comparison = c(rep(comp1$both$Comparison[1], 4), rep(comp2$both$Comparison[1], 4),
                 rep(comp4$both$Comparison[1], 4), rep(comp3$both$Comparison[1], 4))
)

otu_overlap_env$Type <- factor(otu_overlap_env$Type, 
                               levels = c("NewtBd only", "NewtnoBd only", "Pond only", 
                                          "Bislineata only", "Stream only", "Cinereus only", "Forest only", "Shared"))

otu_sum_plot <- ggplot(otu_overlap_env, aes(y=NumOTUs , x=SampleType, fill=Type, order=Type)) +
  ylab("Number of ASVs") + 
  geom_bar(stat="identity", colour="black") +
  facet_grid(. ~ Comparison, scales = "free") +
  scale_fill_manual(limits = levels(otu_overlap_env$Type),
                    values = c("tan4", "blue", "red", "yellow", "green", "deepskyblue3", "gray39", "gray39"))+
  geom_text(aes(label=NumOTUs), position="stack", vjust=-0.3, size=2.8) +
  theme_bw()

ggsave("fullcommunity_env_richnes_asv.pdf", otu_sum_plot, width = 160, height = 90, units="mm", dpi = 300)
