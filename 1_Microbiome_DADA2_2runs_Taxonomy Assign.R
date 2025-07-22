library(dada2)
library(phyloseq)
library(seqRFLP)

# --------- 1. Define Paths to Runs ---------
path.run1 <- "./bacteria_sequences_run1/"
path.run2 <- "./bacteria_sequences_run2/"

list.files(path.run1)
list.files(path.run2)

# --------- 2. Process Run 1 ---------
fnFs1 <- sort(list.files(path.run1, pattern = "_R1_001.fastq", full.names = TRUE))
fnRs1 <- sort(list.files(path.run1, pattern = "_R2_001.fastq", full.names = TRUE))
sample.names1 <- sapply(strsplit(basename(fnFs1), "_"), `[`, 1)

# Quality inspection
plotQualityProfile(fnFs1[1:4])
plotQualityProfile(fnRs1[1:4])

# Filter and trim
filtFs1 <- file.path(path.run1, "filtered", paste0(sample.names1, "_F_filt.fastq.gz"))
filtRs1 <- file.path(path.run1, "filtered", paste0(sample.names1, "_R_filt.fastq.gz"))

out.run1 <- filterAndTrim(fnFs1, filtFs1, fnRs1, filtRs1,
                          truncLen = c(270, 180),
                          maxN = 0, maxEE = c(2, 2),
                          trimLeft = 19, trimRight = 23,
                          truncQ = 2, rm.phix = TRUE, compress = TRUE)

# Error learning
errF1 <- learnErrors(filtFs1, multithread = TRUE)
errR1 <- learnErrors(filtRs1, multithread = TRUE)
plotErrors(errF1, nominalQ = TRUE)

# Dereplication and inference
derepFs1 <- derepFastq(filtFs1, verbose = TRUE)
derepRs1 <- derepFastq(filtRs1, verbose = TRUE)
names(derepFs1) <- sample.names1
names(derepRs1) <- sample.names1

dadaFs1 <- dada(derepFs1, err = errF1, multithread = TRUE)
dadaRs1 <- dada(derepRs1, err = errR1, multithread = TRUE)

mergers1 <- mergePairs(dadaFs1, derepFs1, dadaRs1, derepRs1, verbose = TRUE)
seqtab1 <- makeSequenceTable(mergers1)
saveRDS(seqtab1, "seqtab_run1.rds")

# --------- 3. Process Run 2 ---------
fnFs2 <- sort(list.files(path.run2, pattern = "_R1_001.fastq", full.names = TRUE))
fnRs2 <- sort(list.files(path.run2, pattern = "_R2_001.fastq", full.names = TRUE))
sample.names2 <- sapply(strsplit(basename(fnFs2), "_"), `[`, 1)

# Quality inspection
plotQualityProfile(fnFs2[1:4])
plotQualityProfile(fnRs2[1:4])

# Filter and trim
filtFs2 <- file.path(path.run2, "filtered", paste0(sample.names2, "_F_filt.fastq.gz"))
filtRs2 <- file.path(path.run2, "filtered", paste0(sample.names2, "_R_filt.fastq.gz"))

out.run2 <- filterAndTrim(fnFs2, filtFs2, fnRs2, filtRs2,
                          truncLen = c(270, 180),
                          maxN = 0, maxEE = c(2, 2),
                          trimLeft = 19, trimRight = 23,
                          truncQ = 2, rm.phix = TRUE, compress = TRUE)

# Error learning
errF2 <- learnErrors(filtFs2, multithread = TRUE)
errR2 <- learnErrors(filtRs2, multithread = TRUE)
plotErrors(errF2, nominalQ = TRUE)

# Dereplication and inference
derepFs2 <- derepFastq(filtFs2, verbose = TRUE)
derepRs2 <- derepFastq(filtRs2, verbose = TRUE)
names(derepFs2) <- sample.names2
names(derepRs2) <- sample.names2

dadaFs2 <- dada(derepFs2, err = errF2, multithread = TRUE)
dadaRs2 <- dada(derepRs2, err = errR2, multithread = TRUE)

mergers2 <- mergePairs(dadaFs2, derepFs2, dadaRs2, derepRs2, verbose = TRUE)
seqtab2 <- makeSequenceTable(mergers2)
saveRDS(seqtab2, "seqtab_run2.rds")

# --------- 4. Merge Runs & Remove Chimeras ---------
st1 <- readRDS("seqtab_run1.rds")
st2 <- readRDS("seqtab_run2.rds")

# Merge sequence tables from two runs
st.all <- mergeSequenceTables(st1, st2, repeats = "sum")
table(nchar(getSequences(st.all)))

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(st.all, method = "consensus", multithread = TRUE, verbose = TRUE)

# --------- 5. Assign Taxonomy ---------
taxa <- assignTaxonomy(seqtab.nochim, "./Assign_Taxonomy_Dada2/rdp_train_set_16.fa.gz", multithread = TRUE)

chunk.size <- 5000
taxonomy.species <- do.call(rbind, lapply(
  split(1:nrow(taxa), sort(1:nrow(taxa) %% ceiling(nrow(taxa)/chunk.size))),
  function(x) addSpecies(taxa[x, ], "./Assign_Taxonomy_Dada2/rdp_species_assignment_16.fa.gz")
))

# --------- 6. Build Phyloseq Object and Rename ASVs ---------
ps <- phyloseq(
  otu_table(seqtab.nochim, taxa_are_rows = FALSE),
  tax_table(taxonomy.species)
)

asv.names <- paste0("ASV", seq(ntaxa(ps)))
seqs <- taxa_names(ps)
names(seqs) <- asv.names
taxa_names(ps) <- asv.names

# --------- 7. Export Outputs ---------
# Convert OTU table to matrix
site_species <- as(otu_table(ps), "matrix")
rownames(site_species) <- sapply(strsplit(rownames(site_species), "f"), `[`, 1)
species_site <- t(site_species)

# Convert taxonomy table to matrix
tax <- as(tax_table(ps), "matrix")

# Export files
write.csv(species_site, "Sal_feature_table_2runs.csv")
write.csv(tax, "Sal_taxonomy_2runs.csv")
write.csv(seqs, "Sal_feature_DNAsequences_2runs.csv")

# Convert CSV to FASTA
seq_data <- read.csv("Sal_feature_DNAsequences_2runs.csv", header = TRUE)
dataframe2fas(seq_data, file = "Sal_feature_DNAsequences_2runs.fasta")
