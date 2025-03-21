### statistically validate that TCGA BRCA samples are homogenous ###
library(vegan)

# run PCSF_FULL.R until you get expr_data_filt
rm(PCA_sample_info, STN_samples)


# remove healthy samples
PCA_data <- expr_data_filt$counts_filt
PCA_data <- PCA_data[, colnames(PCA_data) %in% PCA_sample_info_filt$sample[PCA_sample_info_filt$group != "control"]]
sample_info <- PCA_sample_info_filt[PCA_sample_info_filt$group != "control", ]

# convert expression data frame into CPM normalised + transposed matrix
PCA_data <- cpm(as.matrix(PCA_data), log = T)
PCA_data <- t(PCA_data)

pca <- prcomp(PCA_data, scale. = T, center = T)
pca_data <- pca$x

pca_var <- pca$sdev^2
pca_var_perc <- round(pca_var/sum(pca_var)*100, digits = 2)

pca_data <- as.data.frame(pca_data)

# Merge sample mapping with PCA data
pca_data <- merge(pca_data, sample_info, by.x = "row.names", by.y = "sample")

# Create a custom colour palette for stages
library(RColorBrewer)
groups <- unique(sample_info[ ,2])
num_colors <- length(groups)
colours <- brewer.pal(n = num_colors, name = "Dark2")
names(colours) <- groups

# Create the PCA plot with color mapping
library(ggrepel)
ggplot(pca_data, aes(PC1, PC2, color = group, shape = sample_type)) +
  geom_point() +
  #geom_text_repel(aes(label = row.names(pca_data)), size = 3) +  # opt for point labels
  scale_color_manual(values = colours) +
  theme_bw() +
  labs(x = paste0('PC1: ', pca_var_perc[1], ' %'),
       y = paste0('PC2: ', pca_var_perc[2], ' %')) +
  theme(
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 18),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16))




# Choose a subset of principal components (e.g., the first two) to represent the data
pca_scores <- pca_data[, 2:3]

# Compute the Euclidean distance matrix on these PCA scores
distance_matrix <- dist(pca_scores)

# Create a grouping factor for your samples.
# If all samples are tumour, you may simply assign the same group.
groups <- factor(sample_info$group)

# Run the dispersion analysis
dispersion <- betadisper(distance_matrix, groups)
print(dispersion)

# Perform an ANOVA to test if the dispersion is significantly different (if comparing groups)
anova_dispersion <- anova(dispersion)
print(anova_dispersion)

# perform a permutation test for pairwise differences
permutest_disp <- permutest(dispersion)
print(permutest_disp)



# For pairwise correlation, ensure samples are in columns
# assume tumour_data has samples as rows and genes as columns, so we transpose it.
tumour_cor <- expr_data_filt$counts_filt
tumour_cor <- tumour_cor[, colnames(tumour_cor) %in% PCA_sample_info_filt$sample[PCA_sample_info_filt$group != "control"]]
tumour_cor <- cor(tumour_cor)

# Print the correlation matrix
print(tumour_cor)

# Calculate the mean correlation (excluding self-correlations)
mean_corr <- mean(tumour_cor[upper.tri(tumour_cor)])
print(mean_corr)

# Visualize the correlation matrix with a heatmap
heatmap(tumour_cor,
        main = "Pairwise Correlation among Tumour Samples",
        col = colorRampPalette(c("blue", "white", "red"))(20))

low_corr_indices <- which(tumour_cor < mean_corr, arr.ind = TRUE)

# Remove duplicate pairs and self comparisons (only keep one triangle of the matrix)
low_corr_indices <- low_corr_indices[low_corr_indices[,1] < low_corr_indices[,2], ]

# Create a data frame with sample names and the correlation values
low_corr_pairs <- data.frame(
  Sample1 = rownames(tumour_cor)[low_corr_indices[,1]],
  Sample2 = rownames(tumour_cor)[low_corr_indices[,2]],
  Correlation = tumour_cor[low_corr_indices]
)


low_corr_pairs <- merge(low_corr_pairs, sample_info[, 1:2], by.x = "Sample1", by.y = "sample")
low_corr_pairs <- merge(low_corr_pairs, sample_info[, 1:2], by.x = "Sample2", by.y = "sample")

