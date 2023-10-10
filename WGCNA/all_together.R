library(WGCNA)
library(tidyverse)
library(ggrepel)


# combine stages
ben <- read.csv("rna_seq_data/ben_master_df.csv")
stage_I <- read.csv("rna_seq_data/stage_I_master_df.csv")
stage_II <- read.csv("rna_seq_data/stage_II_master_df.csv")
stage_III <- read.csv("rna_seq_data/stage_III_master_df.csv")
stage_IV <- read.csv("rna_seq_data/stage_IV_master_df.csv")

kylie_data <- merge(ben, stage_I, by = "X")
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


# Merge sample-stage mapping with PCA data
pca_data <- merge(pca_data, stage_info, by.x = "row.names", by.y = "Sample")

# Create a custom color palette for stages
stage_colors <- c("ben" = "blue", "stage_I" = "green", "stage_II" = "red", "stage_III" = "purple", "stage_IV" = "orange")

# Create the PCA plot with color mapping
ggplot(pca_data, aes(PC1, PC2, color = Stage)) +
  geom_point() +
  geom_text_repel(aes(label = row.names(pca_data)), size = 3) +  # Adjust the label size here
  scale_color_manual(values = stage_colors) + # Use the custom color palette
  labs(x = paste0('PC1: ', pca_var_perc[1], ' %'),
       y = paste0('PC2: ', pca_var_perc[2], ' %'))


outliers <- c("GAMuT_23091", "GAMuT_41828",  # BEN
              "GAMuT_IC257") # Stage I


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


# create adjacency matrix
soft_power <- 12 # selected based on number of samples

adjacency <- adjacency(wgcna_data, power = soft_power)







TOM <- TOMsimilarity(adjacency, TOMType = "signed")

diss_TOM <- 1 - TOM









