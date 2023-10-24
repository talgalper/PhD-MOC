library(PCSF)
library(igraph)
library(tidyverse)
library(reshape2)
library(biomaRt)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")


mocCNVData <- read.table("~/Downloads/mocCNVData.txt", sep = "\t", header = T)
allCNVDataIDConverted <- read.table("~/Downloads/allCNVDataIDConverted.txt", sep = "\t", header = T)


mocVariantData <- read.csv("~/Downloads/mocVariantData.csv", header = T, na.strings = ".")
mean_rank_by_gene <- aggregate(Consequence_Rank ~ SYMBOL, data = mocVariantData, FUN = mean)

gene_ensembl <- getBM(attributes = c("external_gene_name", "ensembl_gene_id"), 
                      filters = "external_gene_name", 
                      values = mean_rank_by_gene$SYMBOL, 
                      mart = ensembl)

mean_consequence <- merge(gene_ensembl, mean_rank_by_gene, by.x = "external_gene_name", by.y = "SYMBOL")
mean_consequence <- mean_consequence[, -1]


load("~/Desktop/WGCNA_RData_large/differential_weights(full).RData")
differential_weigthts <- as.matrix(differential_weigthts)

edge_list <- melt(differential_weigthts)
colnames(edge_list) <- c("node_1", "node_2", "weight")

# add a column indicating +ve or -ve for later
edge_list$correlation <- ifelse(edge_list$weight > 0, "positive", "negative")

# remove new column for PCSF
edge_list <- edge_list[, -4]

# any negative values will be >1. closer to 0 means more significant co-expression 
edge_list$weight <- 1 - edge_list$weight


start_time <- Sys.time()
interactome <- construct_interactome(edge_list)
elapsed_time <- Sys.time() - start_time
print(elapsed_time)

terminals <- setNames(as.numeric(mean_consequence$Consequence_Rank), mean_consequence$ensembl_gene_id)

start_time <- Sys.time()
set.seed(1234)
subnet <- PCSF_rand(interactome, terminals, n = 50, r = 0.1, w = 2, b = 1, mu = 0.0005)
elapsed_time <- Sys.time() - start_time
print(elapsed_time)

plot.PCSF(kylie_subnet, node_label_cex = 15)



