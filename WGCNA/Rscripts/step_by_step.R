library(WGCNA)


# combine stages
stage_I <- read.csv("rna_seq_data/stage_I_master_df.csv")
stage_II <- read.csv("rna_seq_data/stage_II_master_df.csv")
stage_III <- read.csv("rna_seq_data/stage_III_master_df.csv")
stage_IV <- read.csv("rna_seq_data/stage_IV_master_df.csv")

kylie_data <- merge(stage_I, stage_II, by = "X")
kylie_data <- merge(kylie_data, stage_III, by = "X")
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



outliers <- c( "GAMuT_IC257",
               "GAMuT_34951", "GAMuT_134032", 
               "GAMuT_OV1723")

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

soft_power <- 14
temp_cor <- cor
cor <- WGCNA::cor


adjacency <- adjacency(wgcna_data, power = soft_power)
TOM <- TOMsimilarity(adjacency)
diss_TOM <- 1 - TOM

gene_tree <- hclust(as.dist(diss_TOM), method = "average")
sizeGrWindow(12, 9)
plot(gene_tree, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity",
     labels = F, hang = 0.04)

dynamic_mods <- cutreeDynamic(dendro = gene_tree, distM = diss_TOM,
                              deepSplit = 2, pamRespectsDendro = F,
                              minClusterSize = 30)

dynamic_colours <- labels2colors(dynamic_mods)
sizeGrWindow(8, 6)
plotDendroAndColors(gene_tree, dynamic_colours, "Dynamic Tree Cut",
                    dendroLabels = F, hang = 0.03,
                    addGuide = T, guideHang = 0.05,
                    main = "Gene dendrogram and module colours")

ME_list <- moduleEigengenes(wgcna_data, colors = dynamic_colours)
MEs <- ME_list$eigengenes
ME_diss <- 1 - cor(MEs)
METree <- hclust(as.dist(ME_diss), method = "average")
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
abline(h = 0.25, col = 'red')

merge <- mergeCloseModules(wgcna_data, dynamic_colours, cutHeight = 0.25, verbose = 3)
merged_colours <- merge$colors
merged_MEs <- merge$newMEs

sizeGrWindow(12, 9)
plotDendroAndColors(gene_tree, cbind(dynamic_colours, merged_colours), 
                    c("DynamicTreeCut","Mergeddynamic"), 
                    dendroLabels = F, hang = 0.03, addGuide = T, guideHang = 0.05)



# compare the this run to the original
original <- bwnet$dendrograms[[1]]
original <- original$height

step_by_step <- gene_tree$height
comparison <- setdiff(unlist(original), unlist(step_by_step))
length(comparison)


