library(biomaRt)
library(tidyverse)
library(WGCNA)
library(edgeR)
library(DESeq2)
library(matrixStats)
library(gridExtra)
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


# load data
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

# combine all tumour samples
all_subtypes <- cbind(LumA_unstranded, LumB_unstranded, Her2_unstranded, Basal_unstranded)

# QC + combines tumour and control samples
all_subtype_counts_filt <- filter_low_expr(tumour_matrix = all_subtypes,
                                           control_matrix = GTEx_ENS)

# normalisation
all_wgcna_data <- as.matrix(all_subtype_counts_filt)
all_wgcna_data <- varianceStabilizingTransformation(all_wgcna_data)
all_wgcna_data <- as.data.frame(all_wgcna_data)
all_wgcna_data <- t(all_wgcna_data)

# clean env
rm(LumA_unstranded, LumB_unstranded, Her2_unstranded, 
   Basal_unstranded, normal_unstranded, all_subtype_counts_filt,
   GTEx_ENS)
collectGarbage()

# choose soft thresholding power
unsigned_sft_data <- pick_power(WGCNA_data = all_wgcna_data,
                                network_type = "unsigned")
signed_sft_data <- pick_power(WGCNA_data = all_wgcna_data,
                              network_type = "signed")

# RESTART R AND LOAD WGCNA ONLY
# identify modules: TOMType = "signed", networkType = "unsigned"
all_subtype_bwnet <- network_modules(WGCNA_data = all_wgcna_data,
                                     Power = 4)
save(all_subtype_bwnet, file = "BRCA/RData/all_GTEx/all_subtype_GTEx_bwnet.RData")


# create tumour and control adj matrix
all_adjacencies <- sep_adj_matrix(WGCNA_data = all_wgcna_data,
                                  tumour_expr_df = all_subtypes,
                                  control_expr_df = normal_unstranded,
                                  power = 10)

# save adj matrix
save(all_adjacencies, file = "../../../../Desktop/WGCNA_BRCA_large_files/all_subtype_adj.RData")

load("../../../../Desktop/WGCNA_BRCA_large_files/all_subtype_adj.RData")












multidata <- multiData(Reference = benign_adj, 
                       Test = disease_adj)


multicolour <- list(Reference = bwnet$colors)








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















