# Load the progress package
library(progress)

# Define the number of weight values to test (including 0)
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
            druggability_w * kylie_pcsf_master$druggability -
            citation_w * citation_norm
          
          # Rank the genes based on the combined score
          kylie_pcsf_master_edit <- kylie_pcsf_master
          kylie_pcsf_master_edit$combined_score <- combined_score
          kylie_pcsf_master_edit <- kylie_pcsf_master_edit[order(-kylie_pcsf_master_edit$combined_score), ]
          
          # Select the top 10 genes
          top_genes <- head(kylie_pcsf_master_edit$external_gene_name, 10)
          
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

write.csv(gene_counts, "string_run/kylie_sensitivity_rank.csv")




