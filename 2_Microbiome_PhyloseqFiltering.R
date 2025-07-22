# Load required packages
library(phyloseq)
library(seqRFLP)
library(Biostrings)
library(decontam)
library(vegan)

# -------------------------
# 1. Phyloseq Integration
# -------------------------

# Import data
featureTab <- read.csv("Sal_feature_table_2runs.csv", row.names = 1)
featureTab <- otu_table(featureTab, taxa_are_rows = TRUE)

taxonomy <- tax_table(as.matrix(read.csv("Sal_taxonomy_2runs.csv", row.names = 1)))
meta_data <- sample_data(read.csv("Salamanders_metadata_final_2020_DEC13.csv", row.names = 1))
sequences <- readDNAStringSet("Sal_feature_DNAsequences_2runs.fasta")

# Create phyloseq object
full_dataset <- merge_phyloseq(featureTab, taxonomy, meta_data, sequences)

# -------------------------
# 2. Cleaning and Filtering
# -------------------------

# Keep only microbiome samples and remove failed control
full_dataset <- subset_samples(full_dataset, DATASET_MICROBIOME == "MICROBIOME")
full_dataset <- subset_samples(full_dataset, Sample != "Extr_blankP2_3")

# Remove singleton ASVs
full_dataset <- filter_taxa(full_dataset, function(x) sum(x > 0) > 1, prune=TRUE)

# Remove unwanted taxa
full_dataset <- subset_taxa(full_dataset, !is.na(Phylum) & Family != "Mitochondria" & Order != "Chloroplast" & Kingdom != "Archaea")

# -------------------------
# 3. Remove Contaminants
# -------------------------

# Identify and remove contaminants
sample_data(full_dataset)$is.neg <- sample_data(full_dataset)$Sample_or_Control == "Control Sample"
contamdf <- isContaminant(full_dataset, method="prevalence", neg="is.neg", threshold=0.1)
contaminants <- rownames(contamdf[contamdf$contaminant == TRUE, ])
full_dataset <- prune_taxa(!taxa_names(full_dataset) %in% contaminants, full_dataset)

# Remove control samples
full_dataset <- subset_samples(full_dataset, Sample_or_Control != "Control Sample")

# Remove low abundance ASVs and low-depth samples
full_dataset <- prune_taxa(taxa_sums(full_dataset) >= 10, full_dataset)
full_dataset <- prune_samples(sample_sums(full_dataset) > 2900, full_dataset)

# -------------------------
# 4. Save Unrarefied Dataset
# -------------------------

save(full_dataset, file = "phyloseq_full_dataset_01prev_unrarefy.RData")
write.csv(sample_data(full_dataset), "metadata_salamanders_unrarefy.csv")
write.csv(otu_table(full_dataset), "asv_table_counts_unrarefy.csv")
write.csv(tax_table(full_dataset), "taxonomy_table_filtered_unrarefy.csv")
dataframe2fas(refseq(full_dataset), file = "sequences_unrarefy.fasta")

# -------------------------
# 5. Rarefaction
# -------------------------

full_datasetRare <- rarefy_even_depth(full_dataset, rngseed = 711, replace = FALSE)
save(full_datasetRare, file = "phyloseq_full_dataset_01prev_rarefy.RData")

write.csv(sample_data(full_datasetRare), "metadata_salamanders_rarefy.csv")
write.csv(otu_table(full_datasetRare), "asv_table_counts_rarefy.csv")
write.csv(tax_table(full_datasetRare), "taxonomy_table_filtered_rarefy.csv")
dataframe2fas(refseq(full_datasetRare), file = "sequences_rarefy.fasta")

# -------------------------
# 6. Alpha Diversity Estimates
# -------------------------

t_otu <- t(as(otu_table(full_datasetRare), "matrix"))
AD2 <- vegan::estimateR(t_otu)
write.csv(AD2, "sp_richness_metrics_rarefy.csv")
