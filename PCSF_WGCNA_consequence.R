library(PCSF)
library(igraph)
library(tidyverse)
library(reshape2)
library(biomaRt)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")


mocCNVData <- read.table("~/Downloads/mocCNVData.txt", sep = "\t", header = T)
allCNVDataIDConverted <- read.table("~/Downloads/allCNVDataIDConverted.txt", sep = "\t", header = T)

allCNsegments <- read.csv("~/Downloads/Source data_All CN segments_with IDs.csv", header = T)


mocVariantData <- read.csv("~/Downloads/mocVariantData.csv", header = T, na.strings = ".")
mean_rank_by_gene <- aggregate(Consequence_Rank ~ SYMBOL, data = mocVariantData, FUN = mean)

gene_ensembl <- getBM(attributes = c("external_gene_name", "ensembl_gene_id"), 
                      filters = "external_gene_name", 
                      values = mean_rank_by_gene$SYMBOL, 
                      mart = ensembl)

mean_consequence <- merge(gene_ensembl, mean_rank_by_gene, by.x = "external_gene_name", by.y = "SYMBOL")
mean_consequence <- mean_consequence[, -1]

edge_list <- melt(differential_weights)
colnames(edge_list) <- c("node_1", "node_2", "score")


# any negative values will be >1. closer to 0 means more significant co-expression 
edge_list$score <- 1 - edge_list$score


interactome <- construct_interactome(edge_list)

terminals <- setNames(as.numeric(mean_consequence$Consequence_Rank), mean_consequence$ensembl_gene_id)

# ubuntu
load("~/Desktop/large_git_files/WGCNA/PCSF_data.RData")

start_time <- Sys.time()
set.seed(1234)
subnet <- PCSF_rand(interactome, terminals, n = 10, r = 0.1, w = 2, b = 1, mu = 0.0005)
elapsed_time <- Sys.time() - start_time
print(elapsed_time)

plot.PCSF(subnet, node_label_cex = 15)



