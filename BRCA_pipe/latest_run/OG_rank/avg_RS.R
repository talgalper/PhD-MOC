### perform rank sensitivity with page_rank

library(tidyverse)
library(progress)
library(biomaRt)
library(data.table)

druggability <- read.csv("../Druggability_analysis/data_general/druggability_scores_annot.csv")
load("../ML/RData/full_fpocket_results.RData")
results_master$Score <- (results_master$Score - min(results_master$Score)) / (max(results_master$Score) - min(results_master$Score))
results_master <- results_master[order(-results_master$`Druggability Score`, -results_master$Score), ]
results_master <- results_master[!duplicated(results_master$uniprot_id), ]
druggability <- subset(druggability, select = c("external_gene_name", "druggability", "CP_score", "highest_score"))
druggability <- unique(druggability)


# PubTator3 counts and get total counts per gene
PubTator <- fread("~/OneDrive - RMIT University/PhD/large_git_files/PubTator3/citaiton_counts_ognsmAnnot.csv") # mac
PubTator$tax_id <- as.character(PubTator$tax_id)
PubTator$entrezgene_id <- as.character(PubTator$entrezgene_id)
PubTator[, combined := ifelse(symbol == "", entrezgene_id, symbol)]
PubTator_counts <- PubTator[
  , .(counts = sum(count)),
  by = combined
]
PubTator_counts <- PubTator_counts[order(-PubTator_counts$counts), ]
colnames(PubTator_counts)[1] <- "symbol"


HHnet_enrich <- read.csv("~/Desktop/HHnetEnrich_pageRank.csv")
colnames(HHnet_enrich)[5] <- "degree_centrality"
PCSF <- read.csv("~/Desktop/PCSF_pageRank.csv")

PCSF <- merge.data.table(as.data.table(PCSF), PubTator_counts, by.x = "external_gene_name", by.y = "symbol", all.x = T)
PCSF <- PCSF[order(-PCSF$degree_centrality), ]
PCSF <- merge(PCSF, druggability, by = "external_gene_name", all.x = T)
PCSF <- PCSF[complete.cases(PCSF), ]

HHnet_enrich <- merge.data.table(as.data.table(HHnet_enrich), PubTator_counts, by.x = "external_gene_name", by.y = "symbol", all.x = T)
HHnet_enrich <- HHnet_enrich[order(-HHnet_enrich$degree), ]
HHnet_enrich <- merge(HHnet_enrich, druggability, by = "external_gene_name", all.x = T)
HHnet_enrich <- HHnet_enrich[complete.cases(HHnet_enrich), ]

# normalise data
PCSF$counts_norm <- log10(PCSF$counts)
PCSF[, c(6:10,17)] <- lapply(PCSF[, c(6:10,17)], function(x) {
  (x - min(x)) / (max(x) - min(x))
})

HHnet_enrich$counts_norm <- log10(HHnet_enrich$counts)
HHnet_enrich[, c(5:9,15)] <- lapply(HHnet_enrich[, c(5:9,15)], function(x) {
  (x - min(x)) / (max(x) - min(x))
})



RS_PCSF_pagRank <- rank_sensitivity(PCSF, ensembl, topx = 10, 
                                   features = list(
                                     betweenness = "betweenness",
                                     centrality = "degree_centrality",
                                     druggability = "highest_score",
                                     eigen_centrality = "eigen_centrality",
                                     closeness = "closeness",
                                     page_rank = "page_rank"), 
                                   step = 0.1)


RS_HHnet_enrich_pageRank <- rank_sensitivity(HHnet_enrich, ensembl, topx = 10, 
                                           features = list(
                                             betweenness = "betweenness",
                                             centrality = "degree_centrality",
                                             druggability = "highest_score",
                                             eigen_centrality = "eigen_centrality",
                                             closeness = "closeness",
                                             page_rank = "page_rank"), 
                                           step = 0.1)


targets <- read.csv("../Druggability_analysis/data_general/target_all_dbs.csv")
targets <- targets[, c(2,4)]
targets <- targets[!duplicated(targets$ensembl_gene_id), ]

tem3 <- RS_HHnet_enrich_pageRank[RS_HHnet_enrich_pageRank$external_gene_name %in% targets$drugBank_target, ]
temp4 <- RS_PCSF_pagRank[RS_PCSF_pagRank$external_gene_name %in% targets$drugBank_target, ]



TWS_HHnet_enrich_pageRank <- target_weight_sensitivity(HHnet_enrich, unique(targets$drugBank_target),
                                                       features = list(
                                                         betweenness = "betweenness",
                                                         centrality = "degree_centrality",
                                                         druggability = "highest_score",
                                                         eigen_centrality = "eigen_centrality",
                                                         closeness = "closeness",
                                                         page_rank = "page_rank"), 
                                                       step = 0.1)

TWS_PCSF_pageRank <- target_weight_sensitivity(PCSF, unique(targets$drugBank_target),
                                               features = list(
                                                 betweenness = "betweenness",
                                                 centrality = "degree_centrality",
                                                 druggability = "highest_score",
                                                 eigen_centrality = "eigen_centrality",
                                                 closeness = "closeness",
                                                 page_rank = "page_rank"), 
                                               step = 0.1)





## Average ranking system
combined_rank_sensitivity <- function(input_data,
                                      # Ensembl object for annotation is optional here;
                                      # if needed you can extend the function to retrieve gene descriptions.
                                      features = list(
                                        betweenness = "betweenness",
                                        centrality = "degree_centrality",
                                        druggability = "highest_score",
                                        eigen_centrality = "eigen_centrality",
                                        closeness = "closeness",
                                        citation = "counts_norm"
                                      ),
                                      # Specify which features should be subtracted (e.g. citation counts)
                                      negative_features = c("citation"),
                                      step = 0.1) {
  # Extract feature names and build corresponding weight column names.
  feature_names <- names(features)
  weight_names <- paste0(feature_names, "_w")
  
  # 1. Dynamically generate all weight combinations for the selected features.
  grid_list <- lapply(seq_along(feature_names), function(i) seq(0, 1, step))
  names(grid_list) <- weight_names
  all_combos <- expand.grid(grid_list)
  
  # 2. Keep only the combinations where the weights sum to 1 (within a tiny tolerance).
  valid_combos <- subset(all_combos, abs(rowSums(all_combos) - 1) < 1e-9)
  
  num_combos <- nrow(valid_combos)
  num_genes <- nrow(input_data)
  gene_names <- input_data$external_gene_name
  
  # 3. Create a matrix to store the rank of every gene for each weight combination.
  rank_matrix <- matrix(NA, nrow = num_genes, ncol = num_combos)
  rownames(rank_matrix) <- gene_names
  colnames(rank_matrix) <- paste0("Combo_", seq_len(num_combos))
  
  # Set up a progress bar (optional).
  pb <- progress_bar$new(
    format = "[:bar] :current/:total (:percent) eta: :eta", 
    total = num_combos
  )
  
  # 4. Loop over each valid weight combination.
  for(i in seq_len(num_combos)) {
    pb$tick()
    current_weights <- valid_combos[i, ]
    
    # Compute combined score dynamically based on the selected features.
    # For each feature, multiply the weight by the corresponding input_data column.
    # Use a negative sign for features specified in negative_features.
    combined_score <- rep(0, num_genes)
    for(feature in feature_names) {
      weight_col <- paste0(feature, "_w")
      sign_factor <- ifelse(feature %in% negative_features, -1, 1)
      combined_score <- combined_score +
        sign_factor * current_weights[[weight_col]] * input_data[[ features[[feature]] ]]
    }
    
    # Rank all genes: genes with higher combined scores receive lower (better) ranks.
    order_indices <- order(combined_score, decreasing = TRUE)
    ranks <- integer(num_genes)
    ranks[order_indices] <- seq_along(combined_score)
    
    # Store these ranks in the rank matrix.
    rank_matrix[, i] <- ranks
  }
  
  # 5. Calculate the average rank and variance for each gene over all weight combinations.
  avg_ranks <- rowMeans(rank_matrix)
  var_ranks <- apply(rank_matrix, 1, var)
  
  # Create a summary data frame.
  result_df <- data.frame(
    gene = gene_names,
    avg_rank = avg_ranks,
    rank_variance = var_ranks,
    stringsAsFactors = FALSE
  )
  
  # Order the result by average rank (lower rank = higher priority).
  result_df <- result_df[order(result_df$avg_rank), ]
  
  return(result_df)
}


RS_HHnet_enrich_new <- combined_rank_sensitivity(HHnet_enrich,
                                                 features = list(
                                                   betweenness = "betweenness",
                                                   centrality = "degree_centrality",
                                                   druggability = "highest_score",
                                                   eigen_centrality = "eigen_centrality",
                                                   closeness = "closeness",
                                                   page_rank = "page_rank"
                                                   ), 
                                                 step = 0.1)

rownames(RS_HHnet_enrich_new) <- NULL
temp2 <- RS_HHnet_enrich_new[RS_HHnet_enrich_new$gene %in% targets$drugBank_target, ]




# RS with ML prediction score
HHnet_enrich <- HHnet_enrich_og
HHnet_enrich <- HHnet_enrich[order(-HHnet_enrich$highest_score, HHnet_enrich$degree), ]
HHnet_enrich <- HHnet_enrich[!duplicated(HHnet_enrich$external_gene_name), ]

feature_matrix <- read.csv("~/Desktop/temp.csv", row.names = 1)
feature_matrix <- feature_matrix[, c(1,105:108)]

gene_ids <- getBM(
  attributes = c("uniprot_gn_id", "external_gene_name"),
  filters = "uniprot_gn_id",
  values = feature_matrix$Protein,
  mart = ensembl)

feature_matrix <- merge(gene_ids, feature_matrix, by.x = "uniprot_gn_id", by.y = "Protein", all = T)
feature_matrix <- feature_matrix[order(-feature_matrix$Prediction_Score_rf), ]
rownames(feature_matrix) <- NULL

HHnet_enrich <- merge.data.table(HHnet_enrich, feature_matrix, by = "external_gene_name", all.x = T)
HHnet_enrich <- HHnet_enrich[order(-HHnet_enrich$Prediction_Score_rf), ]
rownames(HHnet_enrich) <- NULL
HHnet_enrich <- unique(HHnet_enrich, by = c("external_gene_name", "degree_centrality"))
HHnet_enrich <- HHnet_enrich[complete.cases(HHnet_enrich), ]

HHnet_enrich$Prediction_Score_rf <- (HHnet_enrich$Prediction_Score_rf - min(HHnet_enrich$Prediction_Score_rf)) / (max(HHnet_enrich$Prediction_Score_rf) - min(HHnet_enrich$Prediction_Score_rf))


RS_HHnet_enrich_ML <- combined_rank_sensitivity(HHnet_enrich,
                                                features = list(
                                                  betweenness = "betweenness",
                                                  centrality = "degree_centrality",
                                                  druggability = "highest_score",
                                                  eigen_centrality = "eigen_centrality",
                                                  closeness = "closeness",
                                                  page_rank = "page_rank",
                                                  ML_pred = "Prediction_Score_rf"), 
                                                step = 0.1)

rownames(RS_HHnet_enrich_ML) <- NULL






# same function as above but trying to include feature contribution/importance
combined_rank_sensitivity_with_contrib <- function(input_data,
                                                   # Optionally, you can supply known targets if needed
                                                   features = list(
                                                     betweenness = "betweenness",
                                                     centrality = "degree_centrality",
                                                     druggability = "highest_score",
                                                     eigen_centrality = "eigen_centrality",
                                                     closeness = "closeness",
                                                     citation = "counts_norm"
                                                   ),
                                                   # Specify which features should be subtracted (e.g. citation counts)
                                                   negative_features = c("citation"),
                                                   step = 0.1) {
  # Extract feature names and build corresponding weight column names.
  feature_names <- names(features)
  weight_names <- paste0(feature_names, "_w")
  
  # 1. Dynamically generate all weight combinations for the selected features.
  grid_list <- lapply(seq_along(feature_names), function(i) seq(0, 1, step))
  names(grid_list) <- weight_names
  all_combos <- expand.grid(grid_list)
  
  # 2. Keep only the combinations where the weights sum to 1 (within a tiny tolerance).
  valid_combos <- subset(all_combos, abs(rowSums(all_combos) - 1) < 1e-9)
  
  num_combos <- nrow(valid_combos)
  num_genes <- nrow(input_data)
  gene_names <- input_data$external_gene_name
  
  # 3. Create a matrix to store the rank of every gene for each weight combination.
  rank_matrix <- matrix(NA, nrow = num_genes, ncol = num_combos)
  rownames(rank_matrix) <- gene_names
  colnames(rank_matrix) <- paste0("Combo_", seq_len(num_combos))
  
  # Also, initialize a list to store contribution matrices (one per feature).
  contrib_list <- list()
  for(feature in feature_names) {
    contrib_list[[feature]] <- matrix(NA, nrow = num_genes, ncol = num_combos)
    rownames(contrib_list[[feature]]) <- gene_names
    colnames(contrib_list[[feature]]) <- paste0("Combo_", seq_len(num_combos))
  }
  
  # Set up a progress bar (optional).
  pb <- progress_bar$new(
    format = "[:bar] :current/:total (:percent) eta: :eta", 
    total = num_combos
  )
  
  # 4. Loop over each valid weight combination.
  for(i in seq_len(num_combos)) {
    pb$tick()
    current_weights <- valid_combos[i, ]
    
    # Compute combined score and store individual contributions.
    combined_score <- rep(0, num_genes)
    
    for(feature in feature_names) {
      weight_col <- paste0(feature, "_w")
      sign_factor <- ifelse(feature %in% negative_features, -1, 1)
      # Partial contribution for this feature:
      partial_contrib <- sign_factor * current_weights[[weight_col]] * input_data[[ features[[feature]] ]]
      combined_score <- combined_score + partial_contrib
      
      # Store the contribution for this feature and weight combination.
      contrib_list[[feature]][, i] <- partial_contrib
    }
    
    # Rank all genes: genes with higher combined scores receive lower (better) ranks.
    order_indices <- order(combined_score, decreasing = TRUE)
    ranks <- integer(num_genes)
    ranks[order_indices] <- seq_along(combined_score)
    
    # Store these ranks.
    rank_matrix[, i] <- ranks
  }
  
  # 5. Calculate the average rank and variance for each gene over all weight combinations.
  avg_ranks <- rowMeans(rank_matrix)
  var_ranks <- apply(rank_matrix, 1, var)
  
  # 6. For each feature, calculate the average contribution (across combinations) per gene.
  avg_contrib <- lapply(feature_names, function(f) {
    rowMeans(contrib_list[[f]])
  })
  names(avg_contrib) <- feature_names
  
  # 7. Build the final result data frame.
  result_df <- data.frame(
    gene = gene_names,
    avg_rank = avg_ranks,
    rank_variance = var_ranks,
    stringsAsFactors = FALSE
  )
  
  # Add average contribution columns for each feature.
  for(feature in feature_names) {
    col_name <- paste0("avg_contrib_", feature)
    result_df[[col_name]] <- avg_contrib[[feature]]
  }
  
  # Order the results by average rank (lower rank = higher priority).
  result_df <- result_df[order(result_df$avg_rank), ]
  
  return(result_df)
}
