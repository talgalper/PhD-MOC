library(PCSF)
library(igraph)
library(tidyverse)
library(reshape2)
library(biomaRt)
library(RobustRankAggreg)
library(gridExtra)
library(enrichR)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")


#mocCNVData <- read.table("~/Downloads/mocCNVData.txt", sep = "\t", header = T)
#allCNVDataIDConverted <- read.table("~/Downloads/allCNVDataIDConverted.txt", sep = "\t", header = T)
#
#allCNsegments <- read.csv("~/Downloads/Source data_All CN segments_with IDs.csv", header = T)


mocVariantData <- read.csv("PCSF/data/mocVariantData.csv", header = T, na.strings = ".")
mean_rank_by_gene <- aggregate(Consequence_Rank ~ SYMBOL, data = mocVariantData, FUN = mean)

gene_ensembl <- getBM(attributes = c("external_gene_name", "ensembl_gene_id"), 
                      filters = "external_gene_name", 
                      values = mean_rank_by_gene$SYMBOL, 
                      mart = ensembl)

mean_consequence <- merge(gene_ensembl, mean_rank_by_gene, by.x = "external_gene_name", by.y = "SYMBOL")
mean_consequence <- mean_consequence[, -1]


#edge_list <- melt(differential_weights)
#colnames(edge_list) <- c("node_1", "node_2", "score")


# any negative values will be >1. closer to 0 means more significant co-expression 
#edge_list$score <- 1 - edge_list$score


interactome <- construct_interactome(edge_list)

terminals <- setNames(as.numeric(mean_consequence$Consequence_Rank), mean_consequence$ensembl_gene_id)

# ubuntu
load("~/Desktop/large_git_files/WGCNA/PCSF_data.RData")
# mac
load("~/OneDrive - RMIT University/PhD/large_git_files/WGCNA/PCSF_data.RData")


start_time <- Sys.time()
set.seed(1234)
subnet <- PCSF_rand(interactome, terminals, n = 10, r = 0.1, w = 2, b = 1, mu = 0.0005)
elapsed_time <- Sys.time() - start_time
print(elapsed_time)

plot.PCSF(subnet, node_label_cex = 15)




load("PCSF/data/PCSF_subnet(n=10).RData")

# extract cluster data
clust <- clusters(subnet)
df <- data.frame(gene_id = names(clust$membership), cluster = factor(clust$membership))
betweenness <- betweenness(subnet) 
centrality <- degree(subnet) 
df$betweenness <- betweenness[as.character(df$gene_id)]
df$degree_centrality <- centrality[as.character(df$gene_id)]
df$betweenness <- as.integer(df$betweenness)
df$degree_centrality <- as.integer(df$degree_centrality)

rownames(df) <- 1:nrow(df)

df <- df[order(-df$degree_centrality), ]


df <- merge(gene_ensembl, df, by.x = "ensembl_gene_id", by.y = "gene_id")

# convert ensembl to uniprot
ensembl_to_uniprot <- getBM(attributes = c("ensembl_gene_id", "uniprot_gn_id"), 
                      filters = "ensembl_gene_id", 
                      values = df$ensembl_gene_id, 
                      mart = ensembl)

# convert external gene name to uniprot
gn_to_uniprot <- getBM(attributes = c("external_gene_name", "uniprot_gn_id"), 
                     filters = "external_gene_name", 
                     values = df$external_gene_name, 
                     mart = ensembl)

PCSF_master <- merge(ensembl_to_uniprot, df, by = "ensembl_gene_id")
PCSF_master2 <- merge(gn_to_uniprot, df, by = "external_gene_name")


# load Fpocket data
af_drugability <- read.csv("../druggability_results/fpocket_druggability.csv")

# merge PCSF data with AF
PCSF_results <- merge(PCSF_master2, af_drugability, by.x = "uniprot_gn_id", by.y = "uniprot_id")

# save results for citation scoring
write.csv(PCSF_results, "PCSF/PCSF_results.csv")

# run citation scoring tool and read in
citation_scores <- read.csv("citation_scores.csv")

PCSF_results <- merge(PCSF_results, citation_scores, by.x = "external_gene_name", by.y = "gene_id")
PCSF_results <- na.omit(PCSF_results) # for some reason the merge created empty duplicates.



# Robust Rank Agrregation
betweenness <- subset(PCSF_results, select = c("ID", "betweenness"))
betweenness <- distinct(betweenness)
betweenness$betweenness <- (betweenness$betweenness - min(betweenness$betweenness)) / (max(betweenness$betweenness) - min(betweenness$betweenness))
betweenness <- betweenness[order(betweenness$betweenness, decreasing = T), ]
betweenness <- betweenness$ID

centrality <- subset(PCSF_results, select = c("ID", "degree_centrality"))
centrality <- distinct(centrality)
centrality$degree_centrality <- (centrality$degree_centrality - min(centrality$degree_centrality)) / (max(centrality$degree_centrality) - min(centrality$degree_centrality))
centrality <- centrality[order(centrality$degree_centrality, decreasing = T), ]
centrality <- centrality$ID

# log transformation of citation scores
citation <- subset(PCSF_results, select = c("ID", "citation_score"))
citation <- distinct(citation)
citation$citation_score <- log(citation$citation_score + 1)
citation$citation_score <- (citation$citation_score - min(citation$citation_score)) / (max(citation$citation_score) - min(citation$citation_score))
citation <- citation[order(citation$citation_score, decreasing = F), ]
citation <- citation$ID

druggability_data <- subset(PCSF_results, select = c("ID", "druggability", "num_drug_pockets"))
druggability_data <- distinct(druggability_data)
druggability <- druggability_data[order(druggability_data$druggability, decreasing = T), ]
druggability <- druggability$ID

num_drug_pockets <- druggability_data[order(druggability_data$num_drug_pockets, decreasing = T), ]
num_drug_pockets <- num_drug_pockets$ID

# add the pocketminer data as another ranking
pocketminer_data <- read.csv("../pocketminer/results/pocketminer_results.csv")
pocketminer_data <- pocketminer_data[pocketminer_data$ID %in% PCSF_results$ID, ]

#pocketminer_gn_id <- getBM(attributes = c("uniprot_gn_id", "external_gene_name"), 
#                           filters = "uniprot_gn_id", 
#                           values = pocketminer_data$uniprot_id, 
#                           mart = ensembl)

#pocketminer_data <- merge(pocketminer_gn_id, pocketminer_data, by.x = "uniprot_gn_id", by.y = "uniprot_id")
#pocketminer_data[pocketminer_data == ""] <- NA
#pocketminer_data <- na.omit(pocketminer_data)
pocketminer_data <- pocketminer_data[order(pocketminer_data$max_hit, decreasing = T), ]
cryptic_pockets <- pocketminer_data$ID

# reorder for num_hits
pocketminer_data <- pocketminer_data[order(pocketminer_data$num_hits, decreasing = T), ]
num_cryp_pockets <- pocketminer_data$ID

rankings <- list(betweenness, centrality, citation, druggability, num_drug_pockets, cryptic_pockets, num_cryp_pockets)


aggregate_ranks <- aggregateRanks(glist = rankings)

# merge gene names back in
temp <- subset(PCSF_results, select = c("external_gene_name", "ID"))
aggregate_ranks <- merge(temp, aggregate_ranks, by.x = "ID", by.y = "Name")
rm(temp)

aggregate_ranks <- aggregate_ranks[order(aggregate_ranks$Score), ]
rownames(aggregate_ranks) <- NULL

# remove scientific notation
options(scipen=999)


write.csv(aggregate_ranks, "PCSF/aggregate_ranks.csv")




## enrichment analysis
aggregate_ranks <- read.csv("PCSF/aggregate_ranks.csv", row.names = 1)

GO_enrichment <- enrichr(aggregate_ranks$external_gene_name, databases = c("GO_Cellular_Component_2023", 
                                                   "GO_Biological_Process_2023",
                                                   "GO_Molecular_Function_2023",
                                                   "KEGG_2021_Human",
                                                   "Proteomics_Drug_Atlas_2023",
                                                   "IDG_Drug_Targets_2022"))

PDA <- GO_enrichment$Proteomics_Drug_Atlas_2023
IDG <- GO_enrichment$IDG_Drug_Targets_2022

CC <- GO_enrichment$GO_Cellular_Component_2023
CC$Category <- "Cellular Component"
CC$Count <- str_count(CC$Genes, ";") + 1
CC <- CC[order(-CC$Count), ]

BP <- GO_enrichment$GO_Biological_Process_2023
BP$Category <- "Biological Process"
BP$Count <- str_count(BP$Genes, ";") + 1
BP <- BP[order(-BP$Count), ]

MF <- GO_enrichment$GO_Molecular_Function_2023
MF$Category <- "Molecular Function"
MF$Count <- str_count(MF$Genes, ";") + 1
MF <- MF[order(-MF$Count), ]


top_terms <- rbind(CC[1:5, ], BP[1:5, ], MF[1:5, ])
top_terms$Term <- gsub("\\s*\\([^)]+\\)", "", top_terms$Term)

# Create a vertical bar plot with grouped x-axis
ggplot(top_terms, aes(x = reorder(Term, Category, FUN = identity), y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(x = "Terms", y = "Counts", fill = "Category") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.margin = margin(l = 150, r = 10, t = 10, b = 10))



KEGG <- GO_enrichment$KEGG_2021_Human
KEGG$Count <- str_count(KEGG$Genes, ";") + 1
KEGG <- KEGG[order(-KEGG$Count), ]

KEGG <- KEGG[order(-KEGG$Count), ]
top_kegg <- KEGG[1:15, ]

ggplot(top_kegg, aes(x = reorder(Term, -Count, FUN = identity), y = Count)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(x = "Terms", y = "Counts", fill = "Category") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))





# because you removed small structures from pocketminer pipeline, you need to remove them from druggability pipeline.
# this code will temporarily fix this but will eventually need to go back and rerun druggabulity package and removing small structures

# need to run on ubuntu
small_structs <- list.files("~/Desktop/pocketminer/results/small_structures/")

PCSF_results <- read.csv("PCSF/PCSF_results.csv", row.names = 1)

IDs <- c()

for (id in small_structs){
  id <- strsplit(id, "-")
  id <- unlist(id)
  id <- paste0(id[1], "-", id[2])
  IDs <- append(IDs, id)
}

PCSF_results_updated <- PCSF_results[!(PCSF_results$ID %in% IDs), ]


plot_data <- pocketminer_data[pocketminer_data$ID %in% PCSF_results$ID, ]

plot_data <- data.frame(num_cryp = plot_data$num_hits,
                        cryp = plot_data$max_hit,
                        num_drug = PCSF_results_updated$num_drug_pockets,
                        drug = PCSF_results_updated$druggability)

num_pockets <- ggplot(plot_data, aes(x = num_drug, y = num_cryp)) +
  geom_point()

pocket_size <- ggplot(plot_data, aes(x = drug, y = cryp)) +
  geom_point()

grid.arrange(num_pockets, pocket_size, ncol = 2)






betweeness_norm <- (PCSF_results$betweenness - min(PCSF_results$betweenness)) / (max(PCSF_results$betweenness) - min(PCSF_results$betweenness))

centrality_norm <- (PCSF_results$degree_centrality - min(PCSF_results$degree_centrality)) / (max(PCSF_results$degree_centrality) - min(PCSF_results$degree_centrality))

# log transformation
citation_norm <- log(PCSF_results$citation_score + 1)
citation_norm <- (citation_norm - min(citation_norm)) / (max(citation_norm) - min(citation_norm))


# Load the progress package
library(progress)

# Define the number of weight values to test
num_weights <- 11

# Create weight value sequences that sum up to 1
betweeness_w_values <- seq(0, 1, length.out = num_weights)
citation_w_values <- seq(0, 1, length.out = num_weights)
centrality_w_values <- seq(0, 1, length.out = num_weights)
druggability_w_values <- seq(0, 1, length.out = num_weights)

# Calculate the total number of iterations
total_iterations <- num_weights ^ 4

# Create a progress bar
pb <- progress_bar$new(total = total_iterations)

# Create an empty data frame with proper column names
sensitivity_results <- data.frame(
  betweeness_weight = numeric(),
  citation_weight = numeric(),
  centrality_weight = numeric(),
  druggability_weight = numeric(),
  top_genes = character(),
  stringsAsFactors = FALSE
)

# Perform sensitivity test
for (betweeness_w in betweeness_w_values) {
  for (citation_w in citation_w_values) {
    for (centrality_w in centrality_w_values) {
      for (druggability_w in druggability_w_values) {
        # Check if the weights sum up to 1
        if (betweeness_w + citation_w + centrality_w + druggability_w == 1) {
          # Combine scores using current weights
          combined_score <- betweeness_w * betweeness_norm +
            centrality_w * centrality_norm +
            druggability_w * PCSF_results$druggability -
            citation_w * citation_norm
          
          # Rank the genes based on the combined score
          PCSF_results_edit <- PCSF_results
          PCSF_results_edit$combined_score <- combined_score
          PCSF_results_edit <- PCSF_results_edit[order(-PCSF_results_edit$combined_score), ]
          
          # Select the top 10 genes
          top_genes <- head(PCSF_results_edit$external_gene_name, 10)
          
          # Add results to the sensitivity_results data frame
          sensitivity_results <- rbind(sensitivity_results, list(
            betweeness_weight = betweeness_w,
            citation_weight = citation_w,
            centrality_weight = centrality_w,
            druggability_weight = druggability_w,
            top_genes = paste(top_genes, collapse = ', ')
          ))
        }
        
        # Increment the progress bar
        pb$tick()
      }
    }
  }
}


# Split the "top_genes" column into a list of genes
sensitivity_results$top_genes_list <- strsplit(sensitivity_results$top_genes, ', ')

# Count the occurrences of each gene
all_genes <- unlist(sensitivity_results$top_genes_list)
all_genes_unique <- unique(all_genes)
count <- sapply(all_genes_unique, function(g) sum(sapply(sensitivity_results$top_genes_list, function(lst) g %in% lst)))

# Create an empty data frame to store gene counts
gene_counts <- data.frame(
  gene = all_genes_unique,
  count = count,
  stringsAsFactors = FALSE
)

# Sort the gene counts by count
gene_counts <- gene_counts[order(-gene_counts$count), ]



counts_when_betweeness_0 <- sensitivity_results[sensitivity_results$betweeness_weight == 0, ]
counts_when_betweeness_0 <- sapply(all_genes_unique, function(g) sum(sapply(counts_when_betweeness_0$top_genes_list, function(lst) g %in% lst)))

gene_counts <- cbind(gene_counts, counts_when_betweeness_0)

counts_when_centrality_0 <- sensitivity_results[sensitivity_results$centrality_weight == 0, ]
counts_when_centrality_0 <- sapply(all_genes_unique, function(g) sum(sapply(counts_when_centrality_0$top_genes_list, function(lst) g %in% lst)))

gene_counts <- cbind(gene_counts, counts_when_centrality_0)

counts_when_citation_0 <- sensitivity_results[sensitivity_results$citation_weight == 0, ]
counts_when_citation_0 <- sapply(all_genes_unique, function(g) sum(sapply(counts_when_citation_0$top_genes_list, function(lst) g %in% lst)))

gene_counts <- cbind(gene_counts, counts_when_citation_0)

counts_when_druggability_0 <- sensitivity_results[sensitivity_results$druggability_weight == 0, ]
counts_when_druggability_0 <- sapply(all_genes_unique, function(g) sum(sapply(counts_when_druggability_0$top_genes_list, function(lst) g %in% lst)))

gene_counts <- cbind(gene_counts, counts_when_druggability_0)

rownames(gene_counts) <- NULL







