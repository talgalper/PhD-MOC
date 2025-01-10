library(tidyverse)
library(progress)
library(biomaRt)

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

colnames(HHnet)[5] <- "degree_centrality"
colnames(HHnet_enrich)[5] <- "degree_centrality"

# column names must match
rank_sensitivity <- function(input_data, ensembl, topx) {
  
  # 1. Generate all weight combinations from 0..1 in steps of 0.1
  all_combos <- expand.grid(
    betweeness_w         = seq(0, 1, 0.1),
    citation_w           = seq(0, 1, 0.1),
    centrality_w         = seq(0, 1, 0.1),
    druggability_w       = seq(0, 1, 0.1),
    eigen_centrality_w   = seq(0, 1, 0.1),
    closeness_w          = seq(0, 1, 0.1)
  )
  
  # 2. Keep only those combos that sum to 1 (± tiny numerical tolerance)
  valid_combos <- subset(
    all_combos, 
    abs(betweeness_w + citation_w + centrality_w + 
          druggability_w + eigen_centrality_w + closeness_w - 1) < 1e-9
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
    
    # Compute combined score
    combined_score <- 
      w$betweeness_w         * input_data$betweenness +
      w$centrality_w         * input_data$degree_centrality +
      w$druggability_w       * input_data$highest_score +
      w$eigen_centrality_w   * input_data$eigen_centrality +
      w$closeness_w          * input_data$closeness -
      w$citation_w           * input_data$counts_norm   # note the minus sign
    
    # Extract top 10 with partial sorting
    top_idx <- order(combined_score, decreasing = TRUE)[1:topx]
    
    top_genes <- input_data$external_gene_name[top_idx]
    
    # Store the result
    results_list[[i]] <- data.frame(
      betweeness_weight       = w$betweeness_w,
      citation_weight         = w$citation_w,
      centrality_weight       = w$centrality_w,
      druggability_weight     = w$druggability_w,
      eigen_centrality_weight = w$eigen_centrality_w,
      closeness_weight        = w$closeness_w,
      top_genes               = paste(top_genes, collapse = ', '),
      stringsAsFactors        = FALSE
    )
  }
  
  # Combine into one data frame
  sensitivity_results <- do.call(rbind, results_list)
  
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
  
  
  weight_cols <- colnames(sensitivity_results[1:(ncol(sensitivity_results)-2)])
  
  # get counts when each variable is 0
  for (w in weight_cols) {
    # "Other" weights = everything except the target w
    other_weights <- setdiff(weight_cols, w)
    
    # Subset rows where the target weight == 0
    # AND none of the other weights == 0
    # i.e. rowSums(...) checks if each row has any zeros among the other weights
    subset_df <- sensitivity_results[
      sensitivity_results[[w]] == 0 &
        rowSums(sensitivity_results[, other_weights] == 0) == 0,
    ]
    
    # check statement
    # cat("Number of count sets in", w, ":", "\n")
    # print(nrow(subset_df))
    
    # Count how many times each gene appears in subset_df$top_genes_list,
    # in the same order as gene_counts$gene
    counts_col <- sapply(gene_counts$gene, function(gene) {
      sum(sapply(subset_df$top_genes_list, function(lst) gene %in% lst))
    })
    
    # Build a descriptive column name, e.g. "counts_when_betweeness_0"
    new_col_name <- paste0("counts_when_", sub("_weight$", "", w), "_0")
    
    # Attach this to gene_counts
    gene_counts[[new_col_name]] <- counts_col
  }
  rownames(gene_counts) <- NULL
  
  # Suppose you want to retrieve descriptions for these genes via biomaRt
  gene_description <- getBM(
    attributes = c("external_gene_name", "description"), 
    filters = "external_gene_name", 
    values = gene_counts$gene,
    mart = ensembl
  )
  gene_description$description <- gsub("\\s*\\[.*?\\]", "", gene_description$description)
  
  final_gene_counts <- merge(gene_description, gene_counts, by.x = "external_gene_name", by.y = "gene", all.y = T)
  
  # Suppose you want to sort by a certain overall count column
  final_gene_counts <- final_gene_counts[order(-final_gene_counts$count), ]
  
  rownames(final_gene_counts) <- NULL
  
  return(final_gene_counts)
}

RS_PCSF <- rank_sensitivity(PCSF, ensembl, topx = 10)
RS_HHnet <- rank_sensitivity(HHnet, ensembl, topx = 10)
RS_HHnet_enrich <- rank_sensitivity(HHnet_enrich, ensembl, topx = 10)



targets <- read.csv("../Druggability_analysis/data_general/target_all_dbs.csv")
targets <- targets[, c(2,4)]
targets <- targets[!duplicated(targets$ensembl_gene_id), ]

temp <- final_gene_counts[final_gene_counts$external_gene_name %in% targets$drugBank_target, ]
temp <- PCSF_input[PCSF_input$external_gene_name %in% targets$drugBank_target, ]
