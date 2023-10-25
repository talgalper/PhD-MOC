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

# subset string interactions from WGCNA
gene_list <- edge_list$node_1
gene_list <- as.list(levels(gene_list))
gene_list <- unlist(gene_list)
write.table(gene_list, "~/Desktop/gene_list.txt", quote = F, row.names = F, col.names = F)

string_edge_data <- read.table("STRING network default edge.csv", header = T, sep = ",", stringsAsFactors = F)
ppi_list <- subset(string_edge_data, select = c("name", "stringdb..score"))
ppi_list <- ppi_list %>% 
  separate(name, sep = " ", into = c("node_1", "del", "node_2"))
ppi_list <- subset(ppi_list, select = c("node_1", "node_2", "stringdb..score"))
ppi_list$node_1 <- gsub(".*.\\.", "", ppi_list$node_1)
ppi_list$node_2 <- gsub(".*.\\.", "", ppi_list$node_2)

string_node_data <- read.table("STRING network default node.csv", header = T, sep = ",", stringsAsFactors = F)
node_list <- subset(string_node_data, select = c("name", "query.term"))
node_list$name <- gsub(".*.\\.", "", node_list$name)
ppi_list$original_order <- seq_len(nrow(ppi_list))
merged_df <- merge(ppi_list, node_list, by.x = "node_1", by.y = "name", all.x = TRUE)
merged_df <- merge(merged_df, node_list, by.x = "node_2", by.y = "name", all.x = TRUE)
merged_df <- merged_df[order(merged_df$original_order), ]

final_df <- merged_df[, c("query.term.x", "query.term.y")]
colnames(final_df) <- c("node_1", "node_2")

subset_edge_list <- edge_list[edge_list$node_1 %in% final$node_1 & edge_list$node_2 %in% unique_interactions$node_2, ]


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



