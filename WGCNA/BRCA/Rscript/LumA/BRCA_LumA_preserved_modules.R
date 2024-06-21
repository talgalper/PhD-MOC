library(biomaRt)
library(tidyverse)
library(WGCNA)
library(edgeR)
library(DESeq2)
library(matrixStats)
library(gridExtra)
library(grid)
library(doParallel)
library(reshape2)
library(igraph)

nCores = 16
registerDoParallel(cores = nCores)
enableWGCNAThreads(nThreads = nCores)
WGCNAnThreads()


# load in data
load("../BRCA_pipe/RData/TCGA_normal.RData")
load("../BRCA_pipe/RData/LumA/DE_data.RData")

# load GTEx data
GTEx_data <- read.table("../BRCA_pipe/gene_reads_2017-06-05_v8_breast_mammary_tissue.gct", skip = 2)
colnames(GTEx_data) <- GTEx_data[1, ]
GTEx_data <- GTEx_data[-1, -1]
rownames(GTEx_data) <- NULL

# opt for having gene Ensembl IDs instead of gene names as rownames (same as TCGA)
GTEx_ENS <- column_to_rownames(GTEx_data, "Name")
rownames(GTEx_ENS) <- gsub("\\.\\d+", "", rownames(GTEx_ENS))
GTEx_ENS <- GTEx_ENS[ , -1]
rownames <- rownames(GTEx_ENS)
GTEx_ENS <- as.data.frame(sapply(GTEx_ENS, as.numeric))
rownames(GTEx_ENS) <- rownames
rm(rownames, GTEx_data)
GTEx_ENS[] <- lapply(GTEx_ENS, function(x){as.integer(x)})


# QC + combines tumour and control samples
counts_filt <- filter_low_expr(tumour_matrix = LumA_unstranded,
                               control_matrix = normal_unstranded,
                               sep = T)

# normalisation
lumA_norm <- vst_norm(counts_filt$tumour)
control_norm <- vst_norm(counts_filt$control)

lumA_intersect <- lumA_norm[, colnames(lumA_norm) %in% colnames(control_norm)]
control_intersect <- control_norm[, colnames(control_norm) %in% colnames(lumA_norm)]

# clean env
rm(counts_filt, LumA_unstranded, normal_unstranded)













# choose soft thresholding power
lumA_sft_data <- pick_power(WGCNA_data = lumA_intersect,
                              network_type = "signed")
control_sft_data <- pick_power(WGCNA_data = control_intersect,
                               network_type = "signed")

# RESTART R AND LOAD WGCNA ONLY
# identify modules: TOMType = "signed", networkType = "unsigned"
lumA_bwnet <- network_modules(WGCNA_data = all_wgcna_data,
                                     Power = 12)
control_bwnet <- network_modules(WGCNA_data = all_wgcna_data,
                                     Power = 11)


save(all_subtype_bwnet, file = "BRCA/RData/all_TCGA/all_subtype_bwnet.RData")
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















