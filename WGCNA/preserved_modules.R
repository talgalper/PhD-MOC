library(WGCNA)
library(tidyverse)
library(ggrepel)
library(edgeR)
library(CorLevelPlot)
library(plyr)
library(gridExtra)
library(biomaRt)
library(matrixStats)
library(reshape2)


# combine stages
ben <- read.csv("rna_seq_data/ben_master_df.csv")
stage_I <- read.csv("rna_seq_data/stage_I_master_df.csv")
stage_II <- read.csv("rna_seq_data/stage_II_master_df.csv")
stage_III <- read.csv("rna_seq_data/stage_III_master_df.csv")
stage_IV <- read.csv("rna_seq_data/stage_IV_master_df.csv")

kylie_data <- merge(ben, stage_I,)
kylie_data <- merge(kylie_data, stage_II, by = "X")
kylie_data <- merge(kylie_data, stage_III, by = "X")
kylie_data <- merge(kylie_data, stage_IV, by = "X")

colnames(kylie_data)[1] <- "gene_id"
rownames(kylie_data) <- kylie_data$gene_id
kylie_data <- kylie_data[, -1] # remove gene ids for calcs


# create stage ID for each sample
stage_info <- data.frame(
  Sample = colnames(kylie_data),
  Stage = character(length(colnames(kylie_data)))
)

stages <- c("ben", "stage_I", "stage_II", "stage_III", "stage_IV")

for (stage in stages) {
  stage_data <- get(stage)
  stage_info$Stage[colnames(kylie_data) %in% colnames(stage_data)] <- stage
}


## removing 0 Variance and low activity genes

# Identify outliers. Searches for 0 variance genes and missing entries
gsg <- goodSamplesGenes(t(kylie_data))
summary(gsg)
gsg$allOK



# top x% of samples
num_genes <- nrow(kylie_data)
top_threshold <- ceiling(num_genes * 0.50)

# min required samples
num_samples <- ncol(kylie_data)
min_required_samples <- ceiling(num_samples * 0.10)



# Extract the gene names and data columns
gene_names <- rownames(kylie_data)
data_columns <- kylie_data

# Initialise an empty data frame
results_df <- data.frame(gene_id = gene_names)

# Iterate through each data column and perform the test
for (col_id in seq_along(data_columns)) {
  # Get the column name and the data
  col_name <- colnames(data_columns)[col_id]
  col_data <- data_columns[, col_id]
  
  # Sort the data and mark genes within the top_threshold rows as TRUE
  sorted_indices <- order(col_data, decreasing = TRUE)
  top_indices <- sorted_indices[1:top_threshold]
  
  # Create a logical vector for gene presence in top_threshold
  gene_presence <- rep(FALSE, nrow(kylie_data))
  gene_presence[top_indices] <- TRUE
  
  # Add to the results data frame
  results_df <- cbind(results_df, gene_presence)
}

# Rename columns 
colnames(results_df) <- c("gene_id", colnames(data_columns))


# Calculate the number of TRUE values for each gene across columns
gene_counts <- rowSums(results_df[, -1]) # Exclude the first column ("gene_id")
failed_genes <- gene_counts < min_required_samples
print(summary(failed_genes))

wgcna_data <- subset(kylie_data, !(rownames(kylie_data) %in% gene_names[failed_genes]))


# logCPM normalisation
wgcna_data <- cpm(wgcna_data, log = T)
wgcna_data <- as.data.frame(wgcna_data)


## data visualisation

# This method visualises data as cluster tree
htree <- hclust(dist(t(wgcna_data)))
plot(htree)


# PCA
pca <- prcomp(t(wgcna_data))
pca_data <- pca$x
pca_var <- pca$sdev^2

pca_var_perc <- round(pca_var/sum(pca_var)*100, digits = 2)

pca_data <- as.data.frame(pca_data)


# Merge sample-stage mapping with PCA data
pca_data <- merge(pca_data, stage_info, by.x = "row.names", by.y = "Sample")

# Create a custom color palette for stages
stage_colors <- c("ben" = "blue", "stage_I" = "green", "stage_II" = "red", "stage_III" = "purple", "stage_IV" = "orange")

# Create the PCA plot with color mapping
ggplot(pca_data, aes(PC1, PC2, color = Stage)) +
  geom_point() +
  geom_text_repel(aes(label = row.names(pca_data)), size = 3) +  # Adjust the label size here
  scale_color_manual(values = stage_colors) + # Use the custom color palette
  theme_bw() +
  labs(x = paste0('PC1: ', pca_var_perc[1], ' %'),
       y = paste0('PC2: ', pca_var_perc[2], ' %'))


#outliers <- c("GAMuT_23091", "GAMuT_41828",  # BEN
#              "GAMuT_IC257") # Stage I


#kylie_data_subset <- kylie_data[, !colnames(kylie_data) %in% outliers]


# separate disease and benign samples for preserved modules function
disease_samples <- c(colnames(stage_I), colnames(stage_II), colnames(stage_III), colnames(stage_IV))
wgcna_disease <- wgcna_data[, colnames(wgcna_data) %in% disease_samples]
wgcna_disease <- t(wgcna_disease)

wgcna_benign <- wgcna_data[, colnames(wgcna_data) %in% colnames(ben)]
wgcna_benign <- t(wgcna_benign)


wgcna_data <- t(wgcna_data)


# Choose a set of soft-threshold powers
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# Call the network topology analysis function
sft <- pickSoftThreshold(wgcna_data,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)

sft_data <- sft$fitIndices


# Visualise to pick power
a1 <- ggplot(sft_data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()

a2 <- ggplot(sft_data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()

grid.arrange(a1, a2, nrow = 2)


# identify modules. this includes both benign and disease groups.
bwnet <- blockwiseModules(wgcna_data,
                          maxBlockSize = 15000,
                          TOMType = "signed",
                          power = 10,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3)


# Plot the dendrogram
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)


# create adjacency matrix
disease_adj <- adjacency(wgcna_disease, power = 10, type = "signed")
benign_adj <- adjacency(wgcna_benign, power = 10, type = "signed")


multidata <- multiData(Reference = benign_adj, 
                       Test = disease_adj)


multicolour <- list(Reference = bwnet$colors)


start_time <- Sys.time()
preserved_modules <- modulePreservation(multiData = multidata,
                                        multiColor = multicolour,
                                        networkType = "signed",
                                        quickCor = 1,
                                        randomSeed = 1234,
                                        verbose = 3,
                                        nPermutations = 10,
                                        testNetworks = 2,
                                        maxModuleSize = max(table(bwnet$colors)),
                                        calculateClusterCoeff = F)
end_time <- Sys.time()
end_time - start_time


# plot modules with median rank and Zsummary threshold

plot_data <- data.frame(
  cluster = rownames(preserved_modules$preservation$Z$ref.Reference$inColumnsAlsoPresentIn.Test),
  moduleSize = preserved_modules$preservation$observed$ref.Reference$inColumnsAlsoPresentIn.Test$moduleSize,
  medianRank.pres = preserved_modules$preservation$observed$ref.Reference$inColumnsAlsoPresentIn.Test$medianRank.pres,
  Zsummary.pres = preserved_modules$preservation$Z$ref.Reference$inColumnsAlsoPresentIn.Test$Zsummary.pres
)


modColors <- unique(plot_data$cluster) 
plotData <-  plot_data[, c(2:ncol(plot_data), 1)]
plotMods <-  !(modColors %in% c("grey", "gold"))

plot1 <- ggplot(plot_data, aes(x = moduleSize, y = medianRank.pres)) +
  geom_point(aes(fill = factor(cluster)), shape = 21, size = 2.4, colour = modColors) +
  scale_x_log10() +
  labs(x = "Module size", y = "Median Rank") +
  theme_minimal() +
  theme(legend.position = "none") +
  geom_text_repel(aes(label = cluster), position = position_nudge(x = 0.1, y = 0.1), color = "black") +
  geom_hline(yintercept = 8, linetype = "dashed") +
  annotate("text", x = 1.5, y = 7.5, label = "Below", size = 3) +
  scale_fill_manual(values = modColors) 


plot2 <- ggplot(plot_data, aes(x = moduleSize, y = Zsummary.pres)) +
  geom_point(aes(fill = factor(cluster)), shape = 21, size = 2.4, colour = modColors) +
  scale_x_log10() +
  labs(x = "Module size", y = "Z Summary") +
  theme_minimal() +
  theme(legend.position = "none") +
  geom_text_repel(aes(label = cluster), position = position_nudge(x = 0.1, y = 0.1), color = "black") +
  geom_hline(yintercept = 10, linetype = "dashed") +
  annotate("text", x = 1.5, y = 11, label = "Above", size = 3) +
  scale_fill_manual(values = modColors)

# Display both plots side by side
grid.arrange(plot1, plot2, ncol = 2)


# extract genes from non-preserved modules
colours <- as.data.frame(bwnet$colors)
colours <- rownames_to_column(colours)
colnames(colours) <- c("genes", "cluster")

non_preserved_modules <- subset(plot_data, medianRank.pres <= 8 & Zsummary.pres >= 10)
non_preserved_genes <- colours[colours$cluster %in% non_preserved_modules$cluster, ]


ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene_names <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), 
                           filters = "ensembl_gene_id", 
                           values = non_preserved_genes$genes, 
                           mart = ensembl)

non_preserved_genes <- cbind(gene_names, non_preserved_genes)
non_preserved_genes <- non_preserved_genes[, -3]

median_rank <- preserved_modules$preservation$observed$ref.Reference$inColumnsAlsoPresentIn.Test
median_rank <- rownames_to_column(median_rank)
median_rank <- subset(median_rank, select = c("rowname", "medianRank.pres"))
merged <- merge(non_preserved_genes, median_rank, by.x = "cluster", by.y = "rowname")

Zsummary <- preserved_modules$preservation$Z$ref.Reference$inColumnsAlsoPresentIn.Test
Zsummary <- rownames_to_column(Zsummary)
Zsummary <- subset(Zsummary, select = c("rowname", "Zsummary.pres"))

merged <- merge(merged, Zsummary, by.x = "cluster", by.y = "rowname")



## diff_i method

# subset genes from adj using non preserved genes
#gene_indices <- rownames(disease_adj) %in% non_preserved_genes$ensembl_gene_id 
#disease_subset <- disease_adj[gene_indices, gene_indices]
#benign_subset <- benign_adj[gene_indices, gene_indices]

# perform diff_i method
sum_matrix <- disease_adj + benign_adj
normalised_scores <- apply(sum_matrix, 2, max)
normalised_scores <- sum_matrix / normalised_scores
median <- rowMedians(normalised_scores)
differential_weights <- normalised_scores - median




## relate module traits

module_eigengenes <- bwnet$MEs


# clean up sample_info and filter for samples in wgcna_data
sample_info <- read.csv("~/Desktop/final_copy/pcsf_kylie/raw_data/All survival_CN_Aug18.csv")
sample_info <- subset(sample_info, select = c("GAMUT_ID", "Grade", "Stage"))

sample_info$GAMUT_ID <- paste0("GAMuT_", sample_info$GAMUT_ID)

wgcna_samples <- rownames(wgcna_data)

sample_info_filtered <- sample_info %>% 
  filter(GAMUT_ID %in% wgcna_samples)

# create list of stage ids 
stage_ids <- c("I","IA","IB","IC", 
               "II", "IIA", "IIB", "IIC", 
               "III", "IIIA", "IIIB", "IIIC", "IIIc", 
               "IV")

# remove unwanted samples
sample_info_filtered <- sample_info_filtered %>%
  filter(is.na(Stage) | Stage %in% stage_ids) %>%
  mutate(Stage = gsub("[AaBbCc]", "", Stage, ignore.case = TRUE)) 

sample_info_filtered$Stage <- ifelse(sample_info_filtered$Grade == "BEN", "BEN", sample_info_filtered$Stage)

# don't need grade column any more
sample_info_filtered <- sample_info_filtered[, -2]
rownames(sample_info_filtered) <- sample_info_filtered$GAMUT_ID


traits <- sample_info_filtered %>% 
  mutate(BEN.vs.all = ifelse(Stage == 'BEN', 1, 0)) %>% 
  select(3)


sample_info_filtered$Stage <- factor(sample_info_filtered$Stage, levels = c("BEN", "I", "II", "III", "IV"))

sample_traits <- binarizeCategoricalColumns(sample_info_filtered$Stage,
                                            includePairwise = FALSE,
                                            includeLevelVsAll = TRUE,
                                            minCount = 1)

traits <- cbind(traits, sample_traits)


# Define numbers of genes and samples
nSamples <- nrow(wgcna_data)
nGenes <- ncol(wgcna_data)


module_trait_corr <- cor(module_eigengenes, traits, use = 'p')
module_trait_corr_pvals <- corPvalueStudent(module_trait_corr, nSamples)


# visualize module-trait association as a heatmap
heatmap_data <- merge(module_eigengenes, traits, by = 'row.names')

heatmap_data <- heatmap_data %>% 
  column_to_rownames(var = 'Row.names')

names(heatmap_data)

CorLevelPlot(heatmap_data,
             x = names(heatmap_data)[24:28],
             y = names(heatmap_data)[1:23],
             col = c("blue1", "skyblue", "white", "pink", "red"))


# pull out a module
module_gene_mapping <- as.data.frame(bwnet$colors)
module <- module_gene_mapping %>% 
  filter(`bwnet$colors` == 'yellow') %>% 
  rownames()


# Calculate the module membership and the associated p-values
# The module membership/intramodular connectivity is calculated as the correlation of the eigengene and the gene expression profile. 
# This quantifies the similarity of all genes on the array to every module.
module_membership_measure <- cor(module_eigengenes, wgcna_data, use = 'p')
module_membership_measure_pvals <- corPvalueStudent(module_membership_measure, nSamples)
module_membership_measure_pvals <- as.data.frame(t(module_membership_measure_pvals))


# Calculate the gene significance and associated p-values 
gene_signf_corr <- cor(wgcna_data, traits$data.IV.vs.all, use = 'p')
gene_signf_corr_pvals <- corPvalueStudent(gene_signf_corr, nSamples)







