library(biomaRt)
library(tidyverse)
library(WGCNA)
library(edgeR)
library(DESeq2)
library(matrixStats)
library(doParallel)
library(reshape2)
library(igraph)

nCores = 8
registerDoParallel(cores = nCores)
enableWGCNAThreads(nThreads = nCores)
WGCNAnThreads()

# load in data
load("../BRCA_pipe/RData/TCGA_normal.RData")
load("../BRCA_pipe/RData/LumA/DE_data.RData")
load("../BRCA_pipe/RData/LumB/DE_data.RData")
load("../BRCA_pipe/RData/Her2/DE_data.RData")
load("../BRCA_pipe/RData/basal/DE_data.RData")

# combine all tumour samples
all_subtypes <- cbind(LumA_unstranded, LumB_unstranded, Her2_unstranded, Basal_unstranded)

# QC + combines tumour and control samples
all_subtype_counts_filt <- filter_low_expr(tumour_matrix = all_subtypes,
                                           control_matrix = normal_unstranded)

# normalisation
all_wgcna_data <- vst_norm(all_subtype_counts_filt)

# plot PCA
control_info <- data.frame(sample = colnames(normal_unstranded),
                           group = rep("control", ncol(normal_unstranded)))
lumA_info <- data.frame(sample = colnames(LumA_unstranded),
                        group = rep("lumA", ncol(LumA_unstranded)))
lumB_info <- data.frame(sample = colnames(LumB_unstranded),
                        group = rep("lumB", ncol(LumB_unstranded)))
her2_info <- data.frame(sample = colnames(Her2_unstranded),
                        group = rep("Her2", ncol(Her2_unstranded)))
basal_info <- data.frame(sample = colnames(Basal_unstranded),
                        group = rep("basal", ncol(Basal_unstranded)))
sample_info <- rbind(control_info, lumA_info, lumB_info, her2_info, basal_info)

PCA_results <- plot_PCA(expr_data = all_wgcna_data,
                        sample_info = sample_info,
                        plot_tree = T,
                        output_plot_data = T)

# identify outliers
dynamicCut <- cutreeDynamic(PCA_results$htree, distM = dist(all_wgcna_data), method = "tree", deepSplit = 2, pamRespectsDendro = FALSE)
outlierSamples <- which(dynamicCut == 0)
cleanExprData <- all_wgcna_data[-outlierSamples, ]
outlierSamples <- all_wgcna_data[outlierSamples, ]
outlierSamples <- sample_info[sample_info$sample %in% rownames(outlierSamples), ]

sample_info_filt <- sample_info[sample_info$sample %in% rownames(cleanExprData), ]

# re-plot PCA with outliers removed
PCA_results_filt <- plot_PCA(expr_data = cleanExprData,
                             sample_info = sample_info_filt,
                             plot_tree = T,
                             output_plot_data = T)

# choose soft thresholding power
sft_data_unsigned <- pick_power(WGCNA_data = cleanExprData,
                                network_type = "unsigned")
sft_data_signed <- pick_power(WGCNA_data = cleanExprData,
                              network_type = "signed")

# identify modules: TOMType = "signed", networkType = "unsigned"
# split data
control_expr <- cleanExprData[rownames(cleanExprData) %in% colnames(normal_unstranded), ]
tumour_expr <- cleanExprData[!rownames(cleanExprData) %in% colnames(normal_unstranded), ]

# clean env
rm(normal_unstranded, LumA_unstranded, LumB_unstranded, Her2_unstranded, Basal_unstranded,
   control_info, lumA_info, lumB_info, her2_info, basal_info,
   all_subtype_counts_filt, all_subtypes, dynamicCut)
collectGarbage()

# RESTART R AND LOAD WGCNA ONLY
library(WGCNA)
library(doParallel)
nCores = 8
registerDoParallel(cores = nCores)
enableWGCNAThreads(nThreads = nCores)
WGCNAnThreads()

control_bwnet <- network_modules(WGCNA_data = control_expr,
                                     Power = 8)
tumour_bwnet <- network_modules(WGCNA_data = tumour_expr,
                                     Power = 8)

save(control_bwnet, file = "BRCA/RData/all_TCGA/control_bwnet.RData")
save(tumour_bwnet, file = "BRCA/RData/all_TCGA/tumour_bwnet.RData")
load("BRCA/RData/all_TCGA/all_subtype_bwnet.RData")

# create tumour and control adj matrix
all_adjacencies <- sep_adj_matrix(WGCNA_data = all_wgcna_data,
                                  tumour_expr_df = all_subtypes,
                                  control_expr_df = normal_unstranded,
                                  power = 10)

# save adj matrix
save(all_adjacencies, file = "../../../../Desktop/WGCNA_BRCA_large_files/all_subtype_adj.RData")
load("../../../../Desktop/WGCNA_BRCA_large_files/all_subtype_adj.RData")

# find preserved modules
multidata <- multiData(Reference = all_adjacencies$control, 
                       Test = all_adjacencies$tumour)

multicolour <- list(Reference = all_subtype_bwnet$colors)

rm(all_adjacencies)
collectGarbage()

start_time <- Sys.time()
preserved_modules <- modulePreservation(multiData = multidata,
                                        multiColor = multicolour,
                                        networkType = "unsigned",
                                        quickCor = 1,
                                        randomSeed = 1234,
                                        verbose = 3,
                                        nPermutations = 10,
                                        testNetworks = 2,
                                        maxModuleSize = max(table(all_subtype_bwnet$colors)),
                                        calculateClusterCoeff = F)

elapsed_time <- Sys.time() - start_time
cat("Elapsed time: ")
print(elapsed_time)












# function for diff_i method
# https://academic.oup.com/bioinformatics/article/36/9/2821/5711285?login=false
diff_i <- function(tumour_adj, control_adj) {
  sum_matrix <- tumour_adj + control_adj
  max_scores <- apply(sum_matrix, 2, max)
  normalised_scores <- sweep(sum_matrix, 2, max_scores, FUN = "/")  
  median <- median(normalised_scores)
  differential_weights <- normalised_scores - median
  
  return(differential_weights)
}

all_subtype_dif_net <- diff_i(tumour_adj = all_adjacencies$tumour,
                              control_adj = all_adjacencies$control)

# clear some memory 
rm(all_adjacencies)
collectGarbage()

all_subtype_edgeList <- melt(all_subtype_dif_net)
colnames(all_subtype_edgeList) <- c("node_1", "node_2", "weight")

edgeList_filt <- all_subtype_edgeList[all_subtype_edgeList$weight > -0.01 & all_subtype_edgeList$weight < 0.01, ]















