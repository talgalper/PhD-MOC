
library(plyr)
library(biomaRt)
library(dplyr) #need to install tidy packages separately because of biomart thing
library(tidyr)
library(tidyverse)
library(progress)
library(rDGIdb)


PCSF_master <- read.csv("intermediate/PCSF_master_unique.csv", row.names = 1)

pocketminer_data <- read.csv("../pocketminer/results/pocketminer_results_3.0.csv")
PCSF_master <- merge(PCSF_master, pocketminer_data, by = "ID")



betweeness_norm <- (PCSF_master$betweenness - min(PCSF_master$betweenness)) / (max(PCSF_master$betweenness) - min(PCSF_master$betweenness))

centrality_norm <- (PCSF_master$degree_centrality - min(PCSF_master$degree_centrality)) / (max(PCSF_master$degree_centrality) - min(PCSF_master$degree_centrality))

# log transformation
citation_norm <- log(PCSF_master$MeSH_count + 1)
citation_norm <- (citation_norm - min(citation_norm)) / (max(citation_norm) - min(citation_norm))

num_drug_pockets_norm <- (PCSF_master$num_drug_pockets - min(PCSF_master$num_drug_pockets)) / (max(PCSF_master$num_drug_pockets) - min(PCSF_master$num_drug_pockets))

num_cryp_pockets_norm <- (PCSF_master$num_hits - min(PCSF_master$num_hits)) / (max(PCSF_master$num_hits) - min(PCSF_master$num_hits))



# Define the number of weight values (genes in top x) to test (including 0)
# i.e. 0, 0.1, 0.2, 0.3, ..., 1
num_weights <- 11

# Create weight value sequences that sum up to 1
betweeness_w_values <- seq(0, 1, length.out = num_weights)
citation_w_values <- seq(0, 1, length.out = num_weights)
centrality_w_values <- seq(0, 1, length.out = num_weights)
druggability_w_values <- seq(0, 1, length.out = num_weights)
cryptic_pocket_w_values <- seq(0, 1, length.out = num_weights)


# Calculate the total number of iterations a.k.a number of variables
total_iterations <- num_weights ^ 5

# Create a progress bar
pb <- progress_bar$new(total = total_iterations)

# Create an empty data frame with proper column names
sensitivity_results <- data.frame(
  betweeness_weight = numeric(),
  citation_weight = numeric(),
  centrality_weight = numeric(),
  druggability_weight = numeric(),
  cryptic_pocket_weight = numeric(),
  top_genes = character(),
  stringsAsFactors = FALSE
)

# Perform sensitivity test
for (betweeness_w in betweeness_w_values) {
  for (citation_w in citation_w_values) {
    for (centrality_w in centrality_w_values) {
      for (druggability_w in druggability_w_values) {
        for (cryptic_pocket_w in cryptic_pocket_w_values) {
          # Check if the weights sum up to 1
          if (betweeness_w + citation_w + centrality_w + druggability_w + cryptic_pocket_w == 1) {
            # Combine scores using current weights
            combined_score <- (betweeness_w * betweeness_norm) +
              (centrality_w * centrality_norm) +
              (druggability_w * PCSF_master$druggability) +
              (cryptic_pocket_w * PCSF_master$max_hit) -
              (citation_w * citation_norm)
            
            # Rank the genes based on the combined score
            PCSF_master_ranked <- PCSF_master
            PCSF_master_ranked$combined_score <- combined_score
            PCSF_master_ranked <- PCSF_master_ranked[order(-PCSF_master_ranked$combined_score), ]
            
            # Select the top 10 genes
            top_genes <- head(PCSF_master_ranked$external_gene_name, 10)
            
            # Add results to the sensitivity_results data frame
            sensitivity_results <- rbind(sensitivity_results, list(
              betweeness_weight = betweeness_w,
              citation_weight = citation_w,
              centrality_weight = centrality_w,
              druggability_weight = druggability_w,
              cryptic_pocket_weight = cryptic_pocket_w,
              top_genes = paste(top_genes, collapse = ', ')
            ))
          }
          
          # Increment the progress bar
          pb$tick()
        }
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


# add scores when variable = 0
counts_when_betweeness_0 <- sensitivity_results[
  sensitivity_results$betweeness_weight == 0 &
    rowSums(sensitivity_results[, c("citation_weight", "centrality_weight", "druggability_weight", "cryptic_pocket_weight")]) != 0, ]

counts_when_betweeness_0 <- sapply(all_genes_unique, function(g) sum(sapply(counts_when_betweeness_0$top_genes_list, function(lst) g %in% lst)))
gene_counts <- cbind(gene_counts, counts_when_betweeness_0)


counts_when_centrality_0 <- sensitivity_results[
  sensitivity_results$centrality_weight == 0 &
    rowSums(sensitivity_results[, c("citation_weight", "betweeness_weight", "druggability_weight", "cryptic_pocket_weight")]) != 0, ]

counts_when_centrality_0 <- sapply(all_genes_unique, function(g) sum(sapply(counts_when_centrality_0$top_genes_list, function(lst) g %in% lst)))
gene_counts <- cbind(gene_counts, counts_when_centrality_0)


counts_when_citation_0 <- sensitivity_results[
  sensitivity_results$citation_weight == 0 &
    rowSums(sensitivity_results[, c("betweeness_weight", "centrality_weight", "druggability_weight", "cryptic_pocket_weight")]) != 0, ]

counts_when_citation_0 <- sapply(all_genes_unique, function(g) sum(sapply(counts_when_citation_0$top_genes_list, function(lst) g %in% lst)))
gene_counts <- cbind(gene_counts, counts_when_citation_0)


counts_when_druggability_0 <- sensitivity_results[
  sensitivity_results$druggability_weight == 0 &
    rowSums(sensitivity_results[, c("citation_weight", "centrality_weight", "betweeness_weight", "cryptic_pocket_weight")]) != 0, ]

counts_when_druggability_0 <- sapply(all_genes_unique, function(g) sum(sapply(counts_when_druggability_0$top_genes_list, function(lst) g %in% lst)))
gene_counts <- cbind(gene_counts, counts_when_druggability_0)


counts_when_cryptic_pockets_0 <- sensitivity_results[
  sensitivity_results$cryptic_pocket_weight == 0 &
    rowSums(sensitivity_results[, c("citation_weight", "centrality_weight", "betweeness_weight", "druggability_weight")]) != 0, ]

counts_when_cryptic_pockets_0 <- sapply(all_genes_unique, function(g) sum(sapply(counts_when_cryptic_pockets_0$top_genes_list, function(lst) g %in% lst)))
gene_counts <- cbind(gene_counts, counts_when_cryptic_pockets_0)



# add gene description
gene_description <- getBM(attributes = c("external_gene_name", "description"), 
                          filters = "external_gene_name", 
                          values = gene_counts$gene, 
                          mart = ensembl)

final_gene_counts <- merge(gene_description, gene_counts, by.x = "external_gene_name", by.y = "gene")

final_gene_counts <- final_gene_counts[order(-final_gene_counts$count), ]

final_gene_counts$description <- gsub("\\s*\\[.*?\\]", "", final_gene_counts$description)


DGIdb <- queryDGIdb(final_gene_counts$external_gene_name)
results <- byGene(DGIdb)
results <- subset(results, select = c("Gene", "DistinctDrugCount"))
final_gene_counts <- merge(final_gene_counts, results, by.x = "external_gene_name", by.y = "Gene")
final_gene_counts <- final_gene_counts[order(-final_gene_counts$count), ]
rownames(final_gene_counts) <- NULL


