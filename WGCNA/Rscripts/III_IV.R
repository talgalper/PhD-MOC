library(WGCNA)
library(edgeR)
library(tidyverse)
library(gridExtra)
library(CorLevelPlot)

# combine stages
benign_data <- read.csv("rna_seq_data/ben_master_df.csv")
stage_III <- read.csv("rna_seq_data/stage_III_master_df.csv")
stage_IV <- read.csv("rna_seq_data/stage_IV_master_df.csv")

kylie_data <- merge(benign_data, stage_III, by = "X")
kylie_data <- merge(kylie_data, stage_IV, by = "X")
colnames(kylie_data)[1] <- "gene_id"
rownames(kylie_data) <- kylie_data$gene_id
#gene_id <- kylie_data$gene_id
kylie_data <- kylie_data[, -1] # remove gene ids for calcs


## Outlier detection

# Identify outliers. Searches for 0 variance genes and missing entries
gsg <- goodSamplesGenes(t(kylie_data))
summary(gsg)
gsg$allOK

# remove any bad genes, if any
bad_genes <- kylie_data[gsg$goodGenes == FALSE, ]
kylie_data <- kylie_data[gsg$goodGenes == TRUE, ]


# This method visualises outliers as cluster tree
htree <- hclust(dist(t(kylie_data)))
plot(htree)

# PCA
pca <- prcomp(t(kylie_data))
pca_data <- pca$x
pca_var <- pca$sdev^2

pca_var_perc <- round(pca_var/sum(pca_var)*100, digits = 2)

pca_data <- as.data.frame(pca_data)

ggplot(pca_data, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca_data)) +
  labs(x = paste0('PC1: ', pca_var_perc[1], ' %'),
       y = paste0('PC2: ', pca_var_perc[2], ' %'))


outliers <- c("GAMuT_41828", "GAMuT_23091")
kylie_data_subset <- kylie_data[, !colnames(kylie_data) %in% outliers]

## Normalisation

# log transformation normalisation
data_log_norm <- log(kylie_data_subset + 1)

# top x% of samples
num_genes <- nrow(data_log_norm)
top_threshold <- ceiling(num_genes * 0.60)

# min required samples
num_samples <- ncol(data_log_norm)
min_required_samples <- ceiling(num_samples * 0.10) # plot sensitivity



# Extract the gene names and data columns
gene_names <- rownames(data_log_norm)
data_columns <- data_log_norm

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
  gene_presence <- rep(FALSE, nrow(data_log_norm))
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

# Remove failed genes from data
wgcna_data <- subset(data_log_norm, !(rownames(data_log_norm) %in% gene_names[failed_genes]))

wgcna_data <- t(wgcna_data)



## Network construction

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

# Convert matrix to numeric
wgcna_data[] <- sapply(wgcna_data, as.numeric)

soft_power <- 18
temp_cor <- cor
cor <- WGCNA::cor


# memory estimate w.r.t blocksize
bwnet <- blockwiseModules(wgcna_data,
                          maxBlockSize = 15000,
                          TOMType = "signed",
                          power = soft_power,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3)


cor <- temp_cor

module_eigengenes <- bwnet$MEs

table(bwnet$colors)

# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)



## relate module traits

# clean up sample_info and filter for samples in wgcna_data
sample_info <- read.csv("~/Desktop/final_copy/pcsf_kylie/raw_data/All survival_CN_Aug18.csv")
sample_info <- subset(sample_info, select = c("GAMUT_ID", "Grade", "Stage"))

sample_info$GAMUT_ID <- paste0("GAMuT_", sample_info$GAMUT_ID)

wgcna_samples <- rownames(wgcna_data)

sample_info_filtered <- sample_info %>% 
  filter(GAMUT_ID %in% wgcna_samples)



# create list of stage ids 
stage_ids <- c("III", "IIIA", "IIIB", "IIIC", "IIIc", 
               "IV")

# remove unwanted samples
sample_info_filtered <- sample_info_filtered %>%
  filter(is.na(Stage) | Stage %in% stage_ids) %>%
  mutate(Stage = gsub("[AaBbCc]", "", Stage, ignore.case = TRUE)) 

# put BEN in stage column where there are NA values
sample_info_filtered <- sample_info_filtered %>%
  mutate(Stage = ifelse(Grade == "BEN", "BEN", Stage))

# don't need grade column any more
sample_info_filtered <- sample_info_filtered[, -2]
rownames(sample_info_filtered) <- sample_info_filtered$GAMUT_ID


traits <- sample_info_filtered %>% 
  mutate(disease_state_bin = ifelse(grepl('BEN', Stage), 0, 1)) %>% 
  select(3)


sample_info_filtered$Stage <- factor(sample_info_filtered$Stage, levels = c("BEN", "III", "IV"))


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

CorLevelPlot(heatmap_data,
             x = names(heatmap_data)[14:16],
             y = names(heatmap_data)[1:13],
             col = c("blue1", "skyblue", "white", "pink", "red"))



module_gene_mapping <- as.data.frame(bwnet$colors)
module_gene_mapping %>% 
  filter(`bwnet$colors` == 'pink') %>% 
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

colours <- as.data.frame(bwnet$colors)
colours <- data.frame(genes = rownames(colours),
                      module = colours$`bwnet$colors`)
colours <- list(Reference = colours) 

data <- adjacency(wgcna_data, power = power, type = "signed")
data <- multiData(Reference = data)

preserved_modules <- modulePreservation(multiData = data,
                                        multiColor = colours,
                                        networkType = "signed")



save(bwnet, wgcna_data, file = "III_IV/WGCNA_III_IV.RData")

