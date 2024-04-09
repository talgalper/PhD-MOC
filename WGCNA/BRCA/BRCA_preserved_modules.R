library(biomaRt)
library(tidyverse)
library(WGCNA)
library(edgeR)


load("../BRCA_pipe/RData/LumA/DE_data.RData")

data <- merge(LumA_unstranded, normal_unstranded, by = "row.names")
data <- data[, -1]

# remove outlier genes
gsg <- goodSamplesGenes(t(data))
summary(gsg)
gsg$allOK

table(gsg$goodGenes)
table(gsg$goodSamples)


data <- data[gsg$goodGenes == TRUE, ]


# plot samples for outlier detection
htree <- hclust(dist(t(data)))
sizeGrWindow(8, 6)
plot(htree, xlab = "", sub = "")

pca <- prcomp(t(data))
pca_data <- pca$x
pca_var <- pca$sdev^2

pca_var_perc <- round(pca_var/sum(pca_var)*100, digits = 2)

pca_data <- as.data.frame(pca_data)

ggplot(pca_data, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca_data)) +
  labs(x = paste0('PC1: ', pca_var_perc[1], ' %'),
       y = paste0('PC2: ', pca_var_perc[2], ' %'))



keepTheseGenes <- (rowSums(cpm(data) > 1) >= ncol(data)/2)
print(summary(keepTheseGenes))

# add gene ids back into df
data <- rownames_to_column(data)

removedGenes <- data$rowname[!keepTheseGenes]
removedGenes <- as.data.frame(removedGenes)
colnames(removedGenes)[1] <- "gene_id"

data <- data[keepTheseGenes, ]
data <- data[, -1]




