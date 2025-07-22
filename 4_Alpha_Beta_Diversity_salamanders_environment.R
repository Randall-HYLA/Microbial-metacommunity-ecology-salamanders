# ===================================
# Load Required Libraries
# ===================================
library(phyloseq)
library(ggplot2)
library(dplyr)
library(stringr)
library(vegan)
library(microbiome)
library(lme4)
library(car)
library(emmeans)

# ===================================
# Load and Filter Phyloseq Object
# ===================================
load("phyloseq_full_dataset_01prev_rarefy.RData")

# Remove sample DMP28 and zero-sum taxa
sal_filt <- subset_samples(full_datasetRare, Sample != "DMP28") %>%
  prune_taxa(taxa_sums(.) > 0, .)

# Create composite variable: Species + Bd
sal_filt@sam_data$species_new <- paste(sal_filt@sam_data$Species, sal_filt@sam_data$Bd)

# ===================================
# Subset to Species of Interest (edit here for each comparison)
# ===================================
# Example: Viridescens with/without Bd vs Pond
sal_only <- subset_samples(sal_filt, Species == "pond" | species_new %in% c("viridescens yes", "viridescens no"))
sal_only <- prune_taxa(taxa_sums(sal_only) > 0, sal_only)

# ===================================
# Sample Summary
# ===================================
sal_only@sam_data %>%
  group_by(Area, Species) %>%
  summarise(n = n())

df_sal_only <- as(sample_data(sal_only), "data.frame")

# ===================================
# Alpha Diversity Calculation
# ===================================
# Observed
observed <- estimate_richness(sal_only, measures = "Observed") %>%
  mutate(alpha = "Observed",
         Sample = as.factor(gsub("X", "", rownames(.))),
         measure = Observed)

# Shannon
shannon <- estimate_richness(sal_only, measures = "Shannon") %>%
  mutate(alpha = "Shannon",
         Sample = as.factor(rownames(.)),
         measure = Shannon)

# Evenness
evenness_vals <- shannon$Shannon / log(observed$Observed)
Evenness <- data.frame(measure = evenness_vals,
                       Sample = shannon$Sample,
                       alpha = "Evenness")

# Combine into one dataframe
alphadiv <- bind_rows(observed[, c("Sample", "measure", "alpha")],
                      shannon[, c("Sample", "measure", "alpha")],
                      Evenness)

# Add metadata
meta <- data.frame(sample_data(sal_only))
alphadiv <- merge(alphadiv, meta, by = "Sample")

# Summary statistics
alphadiv %>%
  group_by(alpha, Species) %>%
  summarise(mean = mean(measure), sd = sd(measure))

# ===================================
# Alpha Diversity Plots
# ===================================
ggplot(subset(alphadiv, alpha == "Shannon"), aes(x = measure)) +
  geom_histogram()

ggplot(subset(alphadiv, alpha == "Observed"), aes(x = measure, y = species_new)) +
  geom_boxplot()

ggplot(subset(alphadiv, alpha == "Evenness"), aes(x = measure, y = species_new)) +
  geom_boxplot()

# ===================================
# Alpha Diversity Modeling
# ===================================
# Observed
mm.observed <- glmer(measure ~ Species + (1 | Area),
                     family = gaussian(link = "log"),
                     data = subset(alphadiv, alpha == "Observed"), na.action = na.omit)
summary(mm.observed)
Anova(mm.observed)
emmeans(mm.observed, pairwise ~ Species)
confint(pairs(emmeans(mm.observed, ~ Species, type = "response")))

# Shannon
mm.shannon <- glmer(measure ~ Species + (1 | Area),
                    family = gaussian(link = "log"),
                    data = subset(alphadiv, alpha == "Shannon"), na.action = na.omit)
summary(mm.shannon)
Anova(mm.shannon)
emmeans(mm.shannon, pairwise ~ Species)

# Evenness
mm.evenness <- glmer(measure ~ Species + (1 | Area),
                     family = gaussian(link = "log"),
                     data = subset(alphadiv, alpha == "Evenness"), na.action = na.omit)
summary(mm.evenness)
Anova(mm.evenness)
emmeans(mm.evenness, pairwise ~ Species)

# ===================================
# Beta Diversity and Dispersion
# ===================================
sals_table <- as.data.frame(otu_table(sal_only))
sals_metadata <- meta(sal_only)
sals_table <- sals_table[, colSums(sals_table) > 0]
sals_table <- as.data.frame(t(sals_table))

# Bray-Curtis Distance
DistBray <- phyloseq::distance(sal_only, method = "bray")

# Dispersion analysis
mod <- betadisper(DistBray, df_sal_only$Species, type = "centroid")
anova(mod)
permutest(mod, permutations = 99, pairwise = TRUE)
mod.HSD <- TukeyHSD(mod)
plot(mod.HSD)

# Distance dataframe
betaBray <- betadisper(DistBray, df_sal_only$Species, type = "median")
dist.Bray <- data.frame(betaBray.distances = betaBray$distances,
                        df_sal_only)

# Summary stats
dist.Bray %>%
  group_by(Species) %>%
  summarise(mean = mean(betaBray.distances),
            sd = sd(betaBray.distances))

# GLM for beta diversity
glm3 <- glm(betaBray.distances ~ Species, data = dist.Bray,
            family = gaussian(link = "log"))
summary(glm3)
Anova(glm3)
emmeans(glm3, pairwise ~ Species)
confint(pairs(emmeans(glm3, ~ Species, type = "response")))
