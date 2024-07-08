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
load("../BRCA_pipe/RData/LumA/DE_data.RData")
load("../BRCA_pipe/RData/LumB/DE_data.RData")
load("../BRCA_pipe/RData/Her2/DE_data.RData")
load("../BRCA_pipe/RData/basal/DE_data.RData")
load("../BRCA_pipe/RData/TCGA_normal.RData")

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

data <- list(GTEx_ENS = GTEx_ENS, 
             normal_unstranded = normal_unstranded, 
             LumA_unstranded = LumA_unstranded, 
             LumB_unstranded = LumB_unstranded, 
             Her2_unstranded = Her2_unstranded, 
             Basal_unstranded = Basal_unstranded)

common_genes <- Reduce(intersect, list(rownames(data$GTEx_ENS), rownames(data$normal_unstranded),
                                       rownames(data$LumA_unstranded), rownames(data$LumB_unstranded), 
                                       rownames(data$Her2_unstranded), rownames(data$Basal_unstranded)))
common_data <- list()
for (i in seq_along(data)) {
  data_common <- data[[i]]
  names <- names(data)
  name <- names[i]
  print(name)
  
  data_common <- data_common[rownames(data_common) %in% common_genes, ]
  common_data[[name]] <- data_common
  
  rm(name, names, data_common, i)
}

all_samples <- cbind(common_data$GTEx_ENS, common_data$normal_unstranded, 
                     common_data$LumA_unstranded, common_data$LumB_unstranded, 
                     common_data$Her2_unstranded, common_data$Basal_unstranded)

tumour_matrix <- cbind( common_data$LumA_unstranded, common_data$LumB_unstranded, 
                        common_data$Her2_unstranded, common_data$Basal_unstranded)

group <- factor(c(rep(1, length(colnames(tumour_matrix))), 
                  rep(2, length(colnames(common_data$normal_unstranded))),
                  rep(3, length(colnames(common_data$GTEx_ENS)))))


counts_filt <- filterByExpr(all_samples, group = group)
print(table(counts_filt))
counts_filt <- all_samples[counts_filt, ]

gsg <- goodSamplesGenes(t(counts_filt))
table(gsg$goodGenes)
table(gsg$goodSamples)

all_data_counts_filt <- vst_norm(counts_filt)



GTEx_info <- data.frame(sample = colnames(common_data$GTEx_ENS),
                        group = rep("GTEx", ncol(common_data$GTEx_ENS)))
TCGA_normal_info <- data.frame(sample = colnames(common_data$normal_unstranded),
                               group = rep("TCGA_normal", ncol(common_data$normal_unstranded)))
lumA_info <- data.frame(sample = colnames(common_data$LumA_unstranded),
                        group = rep("lumA", ncol(common_data$LumA_unstranded)))
lumB_info <- data.frame(sample = colnames(common_data$LumB_unstranded),
                        group = rep("lumB", ncol(common_data$LumB_unstranded)))
her2_info <- data.frame(sample = colnames(common_data$Her2_unstranded),
                        group = rep("Her2", ncol(common_data$Her2_unstranded)))
basal_info <- data.frame(sample = colnames(common_data$Basal_unstranded),
                         group = rep("basal", ncol(common_data$Basal_unstranded)))
sample_info <- rbind(GTEx_info, TCGA_normal_info, lumA_info, lumB_info, her2_info, basal_info)


PCA_results_all <- plot_PCA(expr_data = all_data_counts_filt,
                            sample_info = sample_info,
                            plot_tree = F,
                            output_plot_data = T)

























