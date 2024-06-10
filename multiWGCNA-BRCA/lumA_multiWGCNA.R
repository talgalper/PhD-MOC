library(multiWGCNA)
library(WGCNA)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(tidyverse)
library(DESeq2)
library(edgeR)
library(gridExtra)

# To enable multi-threading
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


# select a thresholding power
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

sft <- pickSoftThreshold(t(wgcna_data),
                         powerVector = power,
                         networkType = "unsigned",
                         verbose = 5)

sft_data <- sft$fitIndices

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


# define variables for `constructNetworks` function
selectedBarcodes <- c(colnames(LumA_unstranded), colnames(normal_unstranded))

sampleTable <- common[ ,-3]
sampleTable <- sampleTable[sampleTable$cases %in% selectedBarcodes, ]

conditions1 = unique(sampleTable[,2])
conditions2 = unique(sampleTable[,3])

startTime <- Sys.time()
# Construct the combined networks and all the sub-networks
LumA_networks <-  constructNetworks(wgcna_data, sampleTable, conditions1, conditions2,
                                    networkType = "unsigned", power = 10,
                                    minModuleSize = 40, maxBlockSize = 25000,
                                    reassignThreshold = 0, minKMEtoStay = 0.7,
                                    mergeCutHeight = 0.10, numericLabels = TRUE,
                                    pamRespectsDendro = FALSE, verbose=3,
                                    saveTOMs = TRUE)

elapsedTime <- Sys.time() - startTime
elapsedTime

# Save results to a list
results=list()
results$overlaps=iterate(autism_networks, overlapComparisons, plot=TRUE)

# Check the reciprocal best matches between the autism and control networks
head(results$overlaps$autism_vs_controls$bestMatches)




