## rank sensitivity without citation score ##


library(tidyverse)
library(progress)
library(biomaRt)

# column names must match
rank_sensitivity <- function(input_data, ensembl, topx) {
  
  # 1. Generate all weight combinations from 0..1 in steps of 0.1 for five features only
  all_combos <- expand.grid(
    betweeness_w = seq(0, 1, 0.1),
    centrality_w = seq(0, 1, 0.1),
    druggability_w = seq(0, 1, 0.1),
    eigen_centrality_w = seq(0, 1, 0.1),
    closeness_w = seq(0, 1, 0.1)
  )
  
  # 2. Keep only those combos that sum to 1 (± tiny numerical tolerance)
  valid_combos <- subset(
    all_combos, 
    abs(betweeness_w + centrality_w + druggability_w + eigen_centrality_w + closeness_w - 1) < 1e-9
  )
  
  # Create a progress bar over valid_combos only
  pb <- progress_bar$new(
    format = "[:bar] :current/:total (:percent) eta: :eta", 
    total = nrow(valid_combos)
  )
  
  # Store results in a list (faster than repeated rbind)
  results_list <- vector("list", length = nrow(valid_combos))
  
  for(i in seq_len(nrow(valid_combos))) {
    pb$tick()
    w <- valid_combos[i, ]
    
    # Compute combined score (without citation penalisation)
    combined_score <- 
      w$betweeness_w * input_data$betweenness +
      w$centrality_w * input_data$degree_centrality +
      w$druggability_w * input_data$highest_score +
      w$eigen_centrality_w * input_data$eigen_centrality +
      w$closeness_w * input_data$closeness
    
    # Extract top x genes
    top_idx <- order(combined_score, decreasing = TRUE)[1:topx]
    top_genes <- input_data$external_gene_name[top_idx]
    
    # Store the result
    results_list[[i]] <- data.frame(
      betweeness_weight = w$betweeness_w,
      centrality_weight = w$centrality_w,
      druggability_weight = w$druggability_w,
      eigen_centrality_weight = w$eigen_centrality_w,
      closeness_weight = w$closeness_w,
      top_genes = paste(top_genes, collapse = ', '),
      stringsAsFactors = FALSE
    )
  }
  
  # Combine into one data frame
  sensitivity_results <- do.call(rbind, results_list)
  sensitivity_results$top_genes_list <- strsplit(sensitivity_results$top_genes, ', ')
  
  # Count the occurrences of each gene
  all_genes <- unlist(sensitivity_results$top_genes_list)
  all_genes_unique <- unique(all_genes)
  count <- sapply(all_genes_unique, function(g) {
    sum(sapply(sensitivity_results$top_genes_list, function(lst) g %in% lst))
  })
  
  gene_counts <- data.frame(
    gene = all_genes_unique,
    count = count,
    stringsAsFactors = FALSE
  )
  
  gene_counts <- gene_counts[order(-gene_counts$count), ]
  
  weight_cols <- colnames(sensitivity_results[1:(ncol(sensitivity_results)-2)])
  
  # Get counts when each variable is 0 (while the other weights are non-zero)
  for (w in weight_cols) {
    other_weights <- setdiff(weight_cols, w)
    subset_df <- sensitivity_results[
      sensitivity_results[[w]] == 0 &
        rowSums(sensitivity_results[, other_weights] == 0) == 0, 
    ]
    
    counts_col <- sapply(gene_counts$gene, function(gene) {
      sum(sapply(subset_df$top_genes_list, function(lst) gene %in% lst))
    })
    
    new_col_name <- paste0("counts_when_", sub("_weight$", "", w), "_0")
    gene_counts[[new_col_name]] <- counts_col
  }
  
  rownames(gene_counts) <- NULL
  
  gene_description <- getBM(
    attributes = c("external_gene_name", "description"), 
    filters = "external_gene_name", 
    values = gene_counts$gene,
    mart = ensembl
  )
  gene_description$description <- gsub("\\s*\\[.*?\\]", "", gene_description$description)
  
  final_gene_counts <- merge(gene_description, gene_counts, by.x = "external_gene_name", by.y = "gene", all.y = TRUE)
  final_gene_counts <- final_gene_counts[order(-final_gene_counts$count), ]
  
  rownames(final_gene_counts) <- NULL
  
  return(final_gene_counts)
}

save(RS_HHnet, RS_HHnet_enrich, RS_PCSF, file = "~/Desktop/temp.RData")


# to run this read in TTD_data from ML script and run part of the TTD_training function
temp <- RS_HHnet_enrich[RS_HHnet_enrich$external_gene_name %in% TTD_clinical_cancer$all_target_genes, ]
temp2 <- RS_HHnet_enrich[RS_HHnet_enrich$external_gene_name %in% TTD_approved_cancer$all_target_genes, ]

temp <- RS_HHnet_enrich_no_cit[RS_HHnet_enrich_no_cit$external_gene_name %in% TTD_clinical_cancer$all_target_genes, ]
temp2 <- RS_HHnet_enrich_no_cit[RS_HHnet_enrich_no_cit$external_gene_name %in% TTD_approved_cancer$all_target_genes, ]

temp <- TTD_clinical_cancer[TTD_clinical_cancer$all_target_genes %in% RS_HHnet_enrich$external_gene_name, ]
temp2 <- TTD_approved_cancer[TTD_approved_cancer$all_target_genes %in% RS_HHnet_enrich$external_gene_name, ]

library(RobustRankAggreg)



RRA <- aggregateRanks()
