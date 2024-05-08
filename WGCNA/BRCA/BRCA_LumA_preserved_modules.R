library(biomaRt)
library(tidyverse)
library(WGCNA)
library(edgeR)
library(DESeq2)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(gridExtra)

# load in data
load("../BRCA_pipe/RData/TCGA_normal.RData")
load("../BRCA_pipe/RData/LumA/lumA_data.RData")
load("../BRCA_pipe/RData/TCGA_query.RData")

data <- merge(LumA_unstranded, normal_unstranded, by = "row.names")
data <- column_to_rownames(data, var = "Row.names")


subtypes <- PanCancerAtlas_subtypes()

common <- merge(clinical_query, subtypes, by.x = "cases", by.y = "pan.samplesID")
common <- subset(common, select = c("cases", "Subtype_Selected", "sample_type", "ajcc_pathologic_stage"))


# remove outlier genes
gsg <- goodSamplesGenes(t(data))
summary(gsg)
gsg$allOK

table(gsg$goodGenes)
table(gsg$goodSamples)


data <- data[gsg$goodGenes == TRUE, ]

# filter low expression genes
group <- factor(c(rep(1, length(colnames(LumA_unstranded))), rep(2, length(colnames(normal_unstranded)))))
counts_filt <- filterByExpr(data, group = group)

# passed genes
table(counts_filt)
data_filt <- data[counts_filt, ]

# normalisation
data_filt <- as.matrix(data_filt)
wgcna_data <- varianceStabilizingTransformation(data_filt)
wgcna_data <- as.data.frame(wgcna_data)



disease_samples <- colnames(LumA_unstranded)
wgcna_disease <- wgcna_data[, colnames(wgcna_data) %in% disease_samples]
wgcna_disease <- t(wgcna_disease)

wgcna_benign <- wgcna_data[, colnames(wgcna_data) %in% colnames(normal_unstranded)]
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


start_time <- Sys.time()
# identify modules. this includes both benign and disease groups.
bwnet <- blockwiseModules(wgcna_data,
                          maxBlockSize = 15000,
                          TOMType = "signed",
                          power = 7,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3)
elapsed_time <- Sys.time() - start_time
system("say run complete")
print(elapsed_time)


# Plot the dendrogram
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)


# create adjacency matrix
disease_adj <- adjacency(wgcna_disease, power = 7, type = "signed")
benign_adj <- adjacency(wgcna_benign, power = 7, type = "signed")


# preserved modules function
multidata <- multiData(Reference = benign_adj, 
                       Test = disease_adj)


multicolour <- list(Reference = bwnet$colors)
system("say run complete")


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
system("say run complete")







# diff_i 
# https://academic.oup.com/bioinformatics/article/36/9/2821/5711285?login=false
sum_matrix <- disease_adj + benign_adj
normalised_scores <- apply(sum_matrix, 2, max)
normalised_scores <- sum_matrix / normalised_scores
median <- rowMedians(normalised_scores)
differential_weights <- normalised_scores - median




## plot samples for outlier detection
#htree <- hclust(dist(t(wgcna_data)), method = "average")
#sizeGrWindow(8, 6)
#plot(htree, xlab = "", sub = "")
#
#
#
## PCA
#pca <- prcomp(t(wgcna_data))
#pca_data <- pca$x
#pca_var <- pca$sdev^2
#
#pca_var_perc <- round(pca_var/sum(pca_var)*100, digits = 2)
#
#pca_data <- as.data.frame(pca_data)
#
#stage_info <- subset(common, select = c("cases", "ajcc_pathologic_stage"))
## Merge sample-stage mapping with PCA data
#pca_data <- merge(pca_data, stage_info, by.x = "row.names", by.y = "cases")
#
## Create a custom color palette for stages
#stage_colors <- c("ben" = "blue", "stage_I" = "green", "stage_II" = "red", "stage_III" = "purple", "stage_IV" = "orange")
#
## Create the PCA plot with color mapping
#ggplot(pca_data, aes(PC1, PC2, color = Stage)) +
#  geom_point() +
#  geom_text_repel(aes(label = row.names(pca_data)), size = 3) +  # Adjust the label size here
#  scale_color_manual(values = stage_colors) + # Use the custom color palette
#  theme_bw() +
#  labs(x = paste0('PC1: ', pca_var_perc[1], ' %'),
#       y = paste0('PC2: ', pca_var_perc[2], ' %'))








