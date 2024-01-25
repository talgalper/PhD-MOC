library(tidyverse)
library(progress)
library(rDGIdb)
library(biomaRt)

data <- read.csv("results/MOC_PCSF_drugability.csv", row.names = 1)
rownames(data) <- NULL
pocketminer_data <- read.csv("../../../pocketminer/results/pocketminer_results_3.0.csv")
master <- merge(data, pocketminer_data, by = "ID")
master <- subset(master, select = c("external_gene_name", "druggability", "struct_score", "max_hit"))



pocketminer_norm <- (master$max_hit - min(master$max_hit)) / 
  (max(master$max_hit) - min(master$max_hit))

structure_norm <- master$struct_score / 100



# Define the number of weight values (genes in top x) to test (including 0)
# i.e. 0, 0.1, 0.2, 0.3, ..., 1
num_weights <- 11

# Create weight value sequences that sum up to 1
fpocket_w_values <- seq(0, 1, length.out = num_weights)
pocketminer_w_values <- seq(0, 1, length.out = num_weights)
structure_w_values <- seq(0, 1, length.out = num_weights)


# Calculate the total number of iterations
total_iterations <- num_weights ^ 3

# Create a progress bar
pb <- progress_bar$new(total = total_iterations)

# Create an empty data frame with proper column names
sensitivity_results <- data.frame(
  fpocket_weight = numeric(),
  pocketminer_weight = numeric(),
  structure_weight = numeric(),
  top_genes = character(),
  stringsAsFactors = FALSE
)

# Perform sensitivity test
for (fpocket_w in fpocket_w_values) {
  for (pocketminer_w in pocketminer_w_values) {
    for (structure_w in structure_w_values) {
        # Check if the weights sum up to 1
        if (fpocket_w + pocketminer_w + structure_w == 1) {
          # Combine scores using current weights
          combined_score <- fpocket_w * master$druggability +
            pocketminer_w * pocketminer_norm +
            structure_w * structure_norm
          
          # Rank the genes based on the combined score
          master_ranked <- master
          master_ranked$combined_score <- combined_score
          master_ranked <- master_ranked[order(-master_ranked$combined_score), ]
          
          # Select the top 10 genes
          top_genes <- head(master_ranked$external_gene_name, 10)
          
          # Add results to the sensitivity_results data frame
          sensitivity_results <- rbind(sensitivity_results, list(
            fpocket_weight = fpocket_w,
            pocketminer_weight = pocketminer_w,
            structure_weight = structure_w,
            top_genes = paste(top_genes, collapse = ', ')
          ))
        }
        
        # Increment the progress bar
        pb$tick()
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



counts_when_fpocket_0 <- sensitivity_results[
  sensitivity_results$fpocket_weight == 0 &
    rowSums(sensitivity_results[, c("fpocket_weight", "pocketminer_weight", "structure_weight")]) != 0, 
]
counts_when_fpocket_0 <- sapply(all_genes_unique, function(g) sum(sapply(counts_when_fpocket_0$top_genes_list, function(lst) g %in% lst)))

gene_counts <- cbind(gene_counts, counts_when_fpocket_0)

counts_when_pocketminer_0 <- sensitivity_results[
  sensitivity_results$pocketminer_weight == 0 &
    rowSums(sensitivity_results[, c("fpocket_weight", "fpocket_weight", "structure_weight")]) != 0, 
]
counts_when_pocketminer_0 <- sapply(all_genes_unique, function(g) sum(sapply(counts_when_pocketminer_0$top_genes_list, function(lst) g %in% lst)))

gene_counts <- cbind(gene_counts, counts_when_pocketminer_0)

counts_when_structure_0 <- sensitivity_results[
  sensitivity_results$structure_weight == 0 &
    rowSums(sensitivity_results[, c("fpocket_weight", "pocketminer_weight", "structure_weight")]) != 0, 
]
counts_when_structure_0 <- sapply(all_genes_unique, function(g) sum(sapply(counts_when_structure_0$top_genes_list, function(lst) g %in% lst)))

final_gene_counts <- cbind(gene_counts, counts_when_structure_0)



# add gene description
gene_description <- getBM(attributes = c("external_gene_name", "description"), 
                          filters = "external_gene_name", 
                          values = final_gene_counts$gene, 
                          mart = ensembl)

final_gene_counts <- merge(gene_description, final_gene_counts, by.x = "external_gene_name", by.y = "gene")

final_gene_counts$description <- gsub("\\s*\\[.*?\\]", "", final_gene_counts$description)


DGIdb <- queryDGIdb(final_gene_counts$gene)
results <- byGene(DGIdb)
results <- subset(results, select = c("Gene", "DistinctDrugCount"))
final_gene_counts <- merge(final_gene_counts, results, by.x = "gene", by.y = "Gene")
final_gene_counts <- final_gene_counts[order(-final_gene_counts$count), ]
rownames(final_gene_counts) <- NULL






cit_scores <- read.csv("intermediate/citation_scores.csv")
master <- merge(master, cit_scores, by.x = "external_gene_name", by.y = "gene_id")

citation_norm <- log(master$citation_score + 1)
citation_norm <- (citation_norm - min(citation_norm)) / (max(citation_norm) - min(citation_norm))

# Create weight value sequences that sum up to 1
fpocket_w_values <- seq(0, 1, length.out = num_weights)
pocketminer_w_values <- seq(0, 1, length.out = num_weights)
structure_w_values <- seq(0, 1, length.out = num_weights)
citation_w_values <- seq(0, 1, length.out = num_weights)


# Calculate the total number of iterations
total_iterations <- num_weights ^ 4

# Create a progress bar
pb <- progress_bar$new(total = total_iterations)

# Create an empty data frame with proper column names
sensitivity_results <- data.frame(
  fpocket_weight = numeric(),
  pocketminer_weight = numeric(),
  structure_weight = numeric(),
  citation_weight = numeric(),
  top_genes = character(),
  stringsAsFactors = FALSE
)

# Perform sensitivity test
for (fpocket_w in fpocket_w_values) {
  for (pocketminer_w in pocketminer_w_values) {
    for (structure_w in structure_w_values) {
      for (citation_w in citation_w_values) {
      # Check if the weights sum up to 1
      if (fpocket_w + pocketminer_w + structure_w + citation_w == 1) {
        # Combine scores using current weights
        combined_score <- fpocket_w * master$druggability +
          pocketminer_w * pocketminer_norm +
          structure_w * structure_norm -
          citation_w * citation_norm
          
        
        # Rank the genes based on the combined score
        master_ranked <- master
        master_ranked$combined_score <- combined_score
        master_ranked <- master_ranked[order(-master_ranked$combined_score), ]
        
        # Select the top 10 genes
        top_genes <- head(master_ranked$external_gene_name, 10)
        
        # Add results to the sensitivity_results data frame
        sensitivity_results <- rbind(sensitivity_results, list(
          fpocket_weight = fpocket_w,
          pocketminer_weight = pocketminer_w,
          structure_weight = structure_w,
          citation_weight = citation_w,
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



counts_when_fpocket_0 <- sensitivity_results[
  sensitivity_results$fpocket_weight == 0 &
    rowSums(sensitivity_results[, c("citation_weight", "pocketminer_weight", "structure_weight")]) != 0, 
]
counts_when_fpocket_0 <- sapply(all_genes_unique, function(g) sum(sapply(counts_when_fpocket_0$top_genes_list, function(lst) g %in% lst)))

gene_counts <- cbind(gene_counts, counts_when_fpocket_0)


counts_when_pocketminer_0 <- sensitivity_results[
  sensitivity_results$pocketminer_weight == 0 &
    rowSums(sensitivity_results[, c("fpocket_weight", "citation_weight", "structure_weight")]) != 0, 
]
counts_when_pocketminer_0 <- sapply(all_genes_unique, function(g) sum(sapply(counts_when_pocketminer_0$top_genes_list, function(lst) g %in% lst)))

gene_counts <- cbind(gene_counts, counts_when_pocketminer_0)


counts_when_structure_0 <- sensitivity_results[
  sensitivity_results$structure_weight == 0 &
    rowSums(sensitivity_results[, c("fpocket_weight", "pocketminer_weight", "citation_weight")]) != 0, 
]
counts_when_structure_0 <- sapply(all_genes_unique, function(g) sum(sapply(counts_when_structure_0$top_genes_list, function(lst) g %in% lst)))

gene_counts <- cbind(gene_counts, counts_when_structure_0)


counts_when_citation_0 <- sensitivity_results[
  sensitivity_results$citation_weight == 0 &
    rowSums(sensitivity_results[, c("fpocket_weight", "pocketminer_weight", "structure_weight")]) != 0, 
]
counts_when_citation_0 <- sapply(all_genes_unique, function(g) sum(sapply(counts_when_citation_0$top_genes_list, function(lst) g %in% lst)))

gene_counts_citation <- cbind(gene_counts, counts_when_citation_0)


# add gene description
gene_description <- getBM(attributes = c("external_gene_name", "description"), 
                          filters = "external_gene_name", 
                          values = gene_counts_citation$gene, 
                          mart = ensembl)

gene_counts_citation <- merge(gene_description, gene_counts_citation, by.x = "external_gene_name", by.y = "gene")

gene_counts_citation$description <- gsub("\\s*\\[.*?\\]", "", gene_counts_citation$description)


DGIdb <- queryDGIdb(gene_counts_citation$gene)
results <- byGene(DGIdb)
results <- subset(results, select = c("Gene", "DistinctDrugCount"))
gene_counts_citation <- merge(gene_counts_citation, results, by.x = "gene", by.y = "Gene")
gene_counts_citation <- gene_counts_citation[order(-gene_counts_citation$count), ]
rownames(gene_counts_citation) <- NULL



save(final_gene_counts, gene_counts_citation, file = "Rscripts/druggability_rank.RData")


