library(multiWGCNA)
library(WGCNA)
library(tidyverse)
library(DESeq2)
library(edgeR)
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
                           group = rep("Control", ncol(GTEx_ENS)))
lumA_info <- data.frame(sample = colnames(LumA_unstranded),
                        group = rep("LumA", ncol(LumA_unstranded)))
lumB_info <- data.frame(sample = colnames(LumB_unstranded),
                        group = rep("LumB", ncol(LumB_unstranded)))
her2_info <- data.frame(sample = colnames(Her2_unstranded),
                        group = rep("Her2", ncol(Her2_unstranded)))
basal_info <- data.frame(sample = colnames(Basal_unstranded),
                         group = rep("Basal", ncol(Basal_unstranded)))
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


wgcna_data_subset <- c(colnames(GTEx_ENS)[1:30],
                       colnames(LumA_unstranded)[1:20],
                       colnames(LumB_unstranded)[1:20],
                       colnames(Her2_unstranded)[1:20],
                       colnames(Basal_unstranded)[1:20])
wgcna_data_subset <- wgcna_data[, colnames(wgcna_data) %in% wgcna_data_subset]
wgcna_data_subset <- wgcna_data_subset[1:2000, 1:ncol(wgcna_data_subset)]
colnames(wgcna_data_subset) <- gsub("-", ".", colnames(wgcna_data_subset)) 


sampleTable <- sample_info
sampleTable$sample <- gsub("-", ".", sampleTable$sample)
sampleTable <- sampleTable[sampleTable$sample %in% colnames(wgcna_data_subset), ]
sampleTable$state <- ifelse(sampleTable$group == "Control", "Healthy", "Tumour")
colnames(sampleTable) <- c("Sample", "Subtype", "State")
sampleTable <- subset(sampleTable, select = c("Sample", "State", "Subtype")) # reorder columns
rownames(sampleTable) <- NULL

conditions1 = sort(unique(sampleTable[,2]))
conditions2 = sort(unique(sampleTable[,3]))

table(sampleTable$Sample %in% colnames(wgcna_data_subset)) # check for miss-matches

# Construct the combined networks and all the sub-networks
multiWGCNA_results <- constructNetworks(wgcna_data_subset, sampleTable, conditions1, conditions2,
                                        networkType = "unsigned", power = 10,
                                        minModuleSize = 40, maxBlockSize = 45000,
                                        reassignThreshold = 0, minKMEtoStay = 0.7,
                                        mergeCutHeight = 0.10, numericLabels = TRUE,
                                        pamRespectsDendro = FALSE, verbose=3,
                                        saveTOMs = FALSE)






