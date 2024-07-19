library(multiWGCNA)
library(WGCNA)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(tidyverse)
library(DESeq2)
library(edgeR)
library(gridExtra)
library(doParallel)

nCores = 8
registerDoParallel(cores = nCores)
enableWGCNAThreads(nThreads = nCores)

# load data
load("../BRCA_pipe/RData/LumA/DE_data.RData")
load("../BRCA_pipe/RData/LumB/DE_data.RData")
load("../BRCA_pipe/RData/Her2/DE_data.RData")
load("../BRCA_pipe/RData/basal/DE_data.RData")
load("../BRCA_pipe/RData/TCGA_normal.RData")

GTEx_data <- read.table("../BRCA_pipe/gene_reads_2017-06-05_v8_breast_mammary_tissue.gct", skip = 2)
colnames(GTEx_data) <- GTEx_data[1, ]
GTEx_data <- GTEx_data[-1, -1]
rownames(GTEx_data) <- NULL

GTEx_ENS <- column_to_rownames(GTEx_data, "Name")
rownames(GTEx_ENS) <- gsub("\\.\\d+", "", rownames(GTEx_ENS))
GTEx_ENS <- GTEx_ENS[ , -1]
rownames <- rownames(GTEx_ENS)
GTEx_ENS <- as.data.frame(sapply(GTEx_ENS, as.numeric))
rownames(GTEx_ENS) <- rownames
rm(rownames, GTEx_data)
GTEx_ENS[] <- lapply(GTEx_ENS, function(x){as.integer(x)})

# subtype sample info
control_info <- data.frame(sample = colnames(GTEx_ENS),
                           group = rep("control", ncol(GTEx_ENS)))
lumA_info <- data.frame(sample = colnames(LumA_unstranded),
                        group = rep("lumA", ncol(LumA_unstranded)))
lumB_info <- data.frame(sample = colnames(LumB_unstranded),
                        group = rep("lumB", ncol(LumB_unstranded)))
her2_info <- data.frame(sample = colnames(Her2_unstranded),
                        group = rep("Her2", ncol(Her2_unstranded)))
basal_info <- data.frame(sample = colnames(Basal_unstranded),
                         group = rep("basal", ncol(Basal_unstranded)))
sample_info <- rbind(control_info, lumA_info, lumB_info, her2_info, basal_info)
rm(lumA_info, lumB_info, her2_info, basal_info, control_info)

# combine all tumour samples
all_subtypes <- cbind(LumA_unstranded, LumB_unstranded, Her2_unstranded, Basal_unstranded)
rm(LumA_unstranded, LumB_unstranded, Her2_unstranded, Basal_unstranded, normal_unstranded)

# filter low expression
counts_filt <- filter_low_expr(tumour_matrix = all_subtypes,
                               control_matrix = GTEx_ENS,
                               sep = T)
# normalise
tumour_norm <- vst_norm(counts_filt$tumour)
control_norm <- vst_norm(counts_filt$control)

# same genes in tumour and common
common_genes <- intersect(colnames(tumour_norm), colnames(control_norm))
tumour_norm <- tumour_norm[, colnames(tumour_norm) %in% common_genes]
control_norm <- control_norm[, colnames(control_norm) %in% common_genes]

wgcna_data <- rbind(control_norm, tumour_norm)
wgcna_data <- as.data.frame(t(wgcna_data))


wgcna_data_subset <- c(colnames(GTEx_ENS)[1:10],
                       colnames(LumA_unstranded)[1:10],
                       colnames(LumB_unstranded)[1:10],
                       colnames(Her2_unstranded)[1:10],
                       colnames(Basal_unstranded)[1:10])
wgcna_data_subset <- wgcna_data[, colnames(wgcna_data) %in% wgcna_data_subset]
wgcna_data_subset <- wgcna_data_subset[1:1000, 1:ncol(wgcna_data_subset)]


sampleTable <- sample_info
sampleTable <- sampleTable[sampleTable$sample %in% colnames(wgcna_data_subset), ]
rownames(sampleTable) <- sampleTable[ , 1]
sampleTable$state <- ifelse(sampleTable$group == "control", "control", "tumour")

colnames(sampleTable) <- c("Sample", "Subtype", "State")
sampleTable <- subset(sampleTable, select = c("Sample", "State", "Subtype"))

sampleTable <- DataFrame(sampleTable)

# define variables for `constructNetworks` function
se <- SummarizedExperiment(assay = list(WGCNA_data = wgcna_data_subset), 
                           colData = list(Sample = sampleTable$Sample,
                                          State = sampleTable$State,
                                          Subtype = sampleTable$Subtype),
                           rowData = rownames(wgcna_data_subset))

conditions1 = unique(sampleTable[,2])
conditions2 = unique(sampleTable[,3])



















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




