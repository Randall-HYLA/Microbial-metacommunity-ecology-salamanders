# Clean and modular R code for microbiome analysis with salamanders

# Load necessary libraries
library(phyloseq)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)

# Cleaned R code for microbiome analysis with salamanders (no custom functions)

# Load phyloseq object
load("phyloseq_full_dataset_01prev_rarefy.RData")

# Filter and clean sample data
sal_filt <- subset_samples(full_datasetRare, Sample != "DMP28")
sal_filt <- prune_taxa(taxa_sums(sal_filt) > 0, sal_filt)
sample_data(sal_filt)$species_new <- paste(sample_data(sal_filt)$Species, sample_data(sal_filt)$Bd)

# --- viridescens Bd+ vs pond ---
vir_bd_pos <- subset_samples(sal_filt, Species == "pond" | species_new == "viridescens yes")
vir_bd_pos <- prune_taxa(taxa_sums(vir_bd_pos) > 0, vir_bd_pos)
vir_bd_pos_asv <- vir_bd_pos %>% tax_glom(taxrank = "ASV") %>% psmelt()
save(vir_bd_pos_asv, file = "viri_bd_pos_stream.RData")

# --- Filter for abundant ASVs ---
vir_bd_pos_filt <- subset(vir_bd_pos_asv, Abundance != 0)
vir_bd_pos_filt$ASV <- as.factor(vir_bd_pos_filt$ASV)
vir_rel_abund <- vir_bd_pos_asv %>%
  group_by(ASV) %>%
  summarize(Abundance = sum(Abundance)) %>%
  mutate(relative_abundance = Abundance / sum(Abundance))

vir_0.001 <- filter(vir_rel_abund, relative_abundance >= 0.001)$ASV

# --- Calculate relative abundance by group ---
asv_grouped <- vir_bd_pos_asv %>%
  select(species_new, ASV, Abundance) %>%
  group_by(species_new, ASV) %>%
  summarize(Abundance = sum(Abundance)) %>%
  filter(Abundance > 0) %>%
  group_by(species_new) %>%
  mutate(relative_abundance = Abundance / sum(Abundance))
asv_filtered <- filter(asv_grouped, ASV %in% vir_0.001)
asv_filtered$duplicate <- duplicated(asv_filtered$ASV)
asv_shared <- filter(asv_filtered, duplicate == TRUE)$ASV
shared_data <- filter(asv_filtered, ASV %in% asv_shared)

# --- Wide format and plot ---
shared_wide <- shared_data %>%
  select(-Abundance, -duplicate) %>%
  pivot_wider(names_from = species_new, values_from = relative_abundance)
ggplot(shared_wide, aes(x = `viridescens yes`, y = `pond NA`)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(x = "N. viridescens Bd+", y = "Pond") +
  xlim(0, 0.01) +
  ylim(0, 0.01) +
  theme_bw()
ggsave("viridescens_bd_corr.pdf", device = "pdf", width = 95, height = 95, units = "mm", dpi = 300)

# --- Correlation ---
cor_data <- shared_data %>%
  select(ASV, species_new, relative_abundance) %>%
  pivot_wider(names_from = species_new, values_from = relative_abundance)
cor.test(cor_data$`pond NA`, cor_data$`viridescens yes`, method = "spearman")

# Repeat for other species comparisons as needed
