# Load Required Libraries
library(phyloseq)
library(ggplot2)
library(dplyr)
library(stringr)
library(vegan)
library(microbiome)
library(lme4)
library(car)
library(sjPlot)
library(emmeans)

# Load Phyloseq Object and Filter Sample
load("phyloseq_full_dataset_01prev_rarefy.RData")
sal_filt <- subset_samples(full_datasetRare, Sample != "DMP28")
sal_filt@sam_data$species_new <- paste(sal_filt@sam_data$Species, sal_filt@sam_data$Bd)
sal_filt <- prune_taxa(taxa_sums(sal_filt) > 0, sal_filt)

# Subset Dataset to Species of Interest
keep_species <- c("viridescens yes", "viridescens no", "cinereus no", "bislineata no")
sal_only <- subset_samples(sal_filt, species_new %in% keep_species)
sal_only <- prune_taxa(taxa_sums(sal_only) > 0, sal_only)

# Prevalence Filtering
source("taxa_summary.R", local = TRUE)
mdt <- fast_melt(sal_only)
prevdt <- mdt[, .(Prevalence = sum(count > 0), TotalCounts = sum(count)), by = TaxaID]
keepTaxa <- prevdt[(Prevalence >= 10 & TotalCounts > 3), TaxaID]
mpra <- transform_sample_counts(sal_only, function(x) x / sum(x))
sal_only <- prune_taxa(keepTaxa, mpra)
sal_only <- prune_taxa(taxa_sums(sal_only) > 0, sal_only)

# Prepare Metadata
df_sal_only <- as(sample_data(sal_only), "data.frame")

# Beta Diversity
sals_table <- t(as.data.frame(otu_table(sal_only)))
sals_table <- sals_table[, colSums(sals_table) > 0]
dist_jaccard <- vegdist(sals_table, method = "jaccard")
dist_bray <- vegdist(sals_table, method = "bray")
mean(dist_bray)

# Ordination Plot
DistBray <- phyloseq::distance(sal_only, method = "jaccard")
ordBray <- ordinate(sal_only, method = "NMDS", distance = DistBray)
plot_ordination(sal_only, ordBray, color = "species_new") +
  geom_point(size = 0.01) +
  stat_ellipse(type = "norm", linetype = 2) +
  geom_text(aes(label = Sample), size = 2, vjust = 1.5) +
  theme_bw()

# PERMDISP Analysis
betaBray <- betadisper(DistBray, df_sal_only$species_new, type = "median")
permutest(betaBray)
anova(betaBray)
boxplot(betaBray$distances ~ betaBray$group)

# GLM on Dispersion Distances
dist_df <- data.frame(distances = betaBray$distances, species_new = df_sal_only$species_new)
glm_model <- glm(distances ~ species_new, data = dist_df, family = gaussian(link = "log"))
summary(glm_model)
Anova(glm_model)

# Mixed Model
mm_model <- glmer(distances ~ species_new + (1 | Area), data = dist_df, family = gaussian(link = "log"))
summary(mm_model)
Anova(mm_model)
plot_model(mm_model, type = "eff", terms = "species_new")

# Alpha Diversity Estimation
observed <- estimate_richness(sal_only, measures = "Observed")
shannon <- estimate_richness(sal_only, measures = "Shannon")
evenness <- shannon$Shannon / log(observed$Observed)

# Combine Alpha Diversity Metrics
alpha_df <- bind_rows(
  observed %>% mutate(alpha = "Observed", measure = Observed),
  shannon %>% mutate(alpha = "Shannon", measure = Shannon),
  data.frame(Sample = rownames(evenness), measure = evenness, alpha = "Evenness")
)

# Merge with Metadata
alpha_df$Sample <- gsub("^X", "", alpha_df$Sample)
alpha_df <- merge(alpha_df, df_sal_only, by = "Sample")

# Plot Evenness
ggplot(subset(alpha_df, alpha == "Evenness"), aes(x = measure, y = species_new)) +
  geom_boxplot()

# Mixed Model on Evenness
mm_evenness <- glmer(measure ~ species_new + (1 | Area), 
                     data = subset(alpha_df, alpha == "Evenness"),
                     family = gaussian(link = "log"))
summary(mm_evenness)
Anova(mm_evenness)

# Save OTU and Taxonomy Tables
write.csv(as.data.frame(otu_table(sal_only)), "Filtered_OTU_Table.csv")
write.csv(as.data.frame(tax_table(sal_only)), "Filtered_Taxonomy.csv")


