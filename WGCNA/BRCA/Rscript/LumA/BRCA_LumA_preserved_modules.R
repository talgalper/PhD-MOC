library(biomaRt)
library(tidyverse)
library(WGCNA)
library(edgeR)
library(DESeq2)
library(matrixStats)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(gridExtra)
library(doParallel)

nCores = 8
registerDoParallel(cores = nCores)
enableWGCNAThreads(nThreads = nCores)

# load in data
load("../BRCA_pipe/RData/TCGA_normal.RData")
load("../BRCA_pipe/RData/LumA/DE_data.RData")
load("../BRCA_pipe/RData/TCGA_query.RData")

# merge normal and disease samples
data <- merge(LumA_unstranded, normal_unstranded, by = "row.names")
data <- column_to_rownames(data, var = "Row.names")

query_output <- getResults(query_TCGA)
clinical_query <- clinical[complete.cases(clinical$ajcc_pathologic_stage), ]
clinical_query <- merge(query_output, clinical_query, by.x = "cases.submitter_id", by.y = "submitter_id")
clinical_query <- subset(clinical_query, select = c("cases", "cases.submitter_id", "ajcc_pathologic_stage", 
                                                    "tissue_or_organ_of_origin", "sample_type"))
# add subtypes to clinical data
common <- merge(clinical_query, subtypes, by.x = "cases", by.y = "pan.samplesID")
common <- subset(common, select = c("cases", "Subtype_Selected", "sample_type", "ajcc_pathologic_stage"))


# remove outlier genes
gsg <- goodSamplesGenes(t(data))
summary(gsg)
gsg$allOK

table(gsg$goodGenes)
table(gsg$goodSamples)

# filter out bad genes. There were no bad samples
data <- data[gsg$goodGenes == TRUE, ]

# filter low expression genes using edgeR function `filterByExpr`
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
                         networkType = "unsigned",
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


## chatGPT reckons pwr 6 instead of 5
start_time <- Sys.time()
# identify modules. this includes both benign and disease groups.
bwnet <- blockwiseModules(wgcna_data,
                          maxBlockSize = 15000,
                          TOMType = "unsigned",
                          power = 6,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3,
                          saveTOMs = FALSE)
elapsed_time <- Sys.time() - start_time
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
#system("say run complete") 


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
#system("say run complete")




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






# diff_i 
# https://academic.oup.com/bioinformatics/article/36/9/2821/5711285?login=false
sum_matrix <- disease_adj + benign_adj
normalised_scores <- apply(sum_matrix, 2, max)
normalised_scores <- sum_matrix / normalised_scores
median <- rowMedians(normalised_scores)
differential_weights <- normalised_scores - median




# plot samples for outlier detection
htree <- hclust(dist(t(data_filt)), method = "average")
sizeGrWindow(8, 6)
plot(htree, xlab = "", sub = "")



# PCA using raw counts
pca <- prcomp(t(data_filt))
pca_data <- pca$x
pca_var <- pca$sdev^2

pca_var_perc <- round(pca_var/sum(pca_var)*100, digits = 2)

pca_data <- as.data.frame(pca_data)

subtype_info <- subset(common, select = c("cases", "Subtype_Selected"))
# Merge sample-stage mapping with PCA data
pca_data <- merge(pca_data, subtype_info, by.x = "row.names", by.y = "cases")

# Create a custom color palette for stages
subtype_colors <- c("BRCA.Normal" = "green", "BRCA.LumA" = "red")

# Create the PCA plot with color mapping
ggplot(pca_data, aes(PC1, PC2, colour = Subtype_Selected)) +
  geom_point() +
  geom_text_repel(aes(label = row.names(pca_data)), size = 3) +  # Adjust the label size here
  scale_color_manual(values = subtype_colors) + # Use the custom color palette
  theme_bw() +
  labs(x = paste0('PC1: ', pca_var_perc[1], ' %'),
       y = paste0('PC2: ', pca_var_perc[2], ' %'))



# PCA using normalised counts (variance stabilised)
pca <- prcomp(wgcna_data)
pca_data <- pca$x
pca_var <- pca$sdev^2

pca_var_perc <- round(pca_var/sum(pca_var)*100, digits = 2)

pca_data <- as.data.frame(pca_data)

subtype_info <- subset(common, select = c("cases", "Subtype_Selected"))
# Merge sample-stage mapping with PCA data
pca_data <- merge(pca_data, subtype_info, by.x = "row.names", by.y = "cases")

# Create a custom color palette for stages
subtype_colors <- c("BRCA.Normal" = "green", "BRCA.LumA" = "red")

# Create the PCA plot with color mapping
ggplot(pca_data, aes(PC1, PC2, colour = Subtype_Selected)) +
  geom_point() +
  geom_text_repel(aes(label = row.names(pca_data)), size = 3) +  # Adjust the label size here
  scale_color_manual(values = subtype_colors) + # Use the custom color palette
  theme_bw() +
  labs(x = paste0('PC1: ', pca_var_perc[1], ' %'),
       y = paste0('PC2: ', pca_var_perc[2], ' %'))





