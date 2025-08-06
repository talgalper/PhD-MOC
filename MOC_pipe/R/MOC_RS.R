### perform rank sensitivity with page_rank ###

library(progress)
library(biomaRt)
library(data.table)
source("R/functions.R")

druggability <- read.csv("../Druggability_analysis/data_general/druggability_scores_annot.csv")
load("../ML/RData/full_fpocket_results.RData")
results_master$Score <- (results_master$Score - min(results_master$Score)) / (max(results_master$Score) - min(results_master$Score))
results_master <- results_master[order(-results_master$`Druggability Score`, -results_master$Score), ]
results_master <- results_master[!duplicated(results_master$uniprot_id), ]
druggability <- subset(druggability, select = c("external_gene_name", "druggability", "CP_score", "highest_score"))
druggability <- unique(druggability)


# PubTator3 counts and get total counts per gene
PubTator <- fread("~/OneDrive - RMIT University/PhD/large_git_files/PubTator3/citaiton_counts_ognsmAnnot.csv")
PubTator$tax_id <- as.character(PubTator$tax_id)
PubTator$entrezgene_id <- as.character(PubTator$entrezgene_id)
PubTator[, combined := ifelse(symbol == "", entrezgene_id, symbol)]
PubTator_counts <- PubTator[
  , .(counts = sum(count)),
  by = combined
]
PubTator_counts <- PubTator_counts[order(-PubTator_counts$counts), ]
colnames(PubTator_counts)[1] <- "symbol"


HHnet_enrich <- read.csv("../Hierarchical_HotNet/MOC/results/df_subnetNeighs.csv")
colnames(HHnet_enrich)[5] <- "degree_centrality"
# PCSF <- read.csv("~/Desktop/PCSF_pageRank.csv")

# PCSF <- merge.data.table(as.data.table(PCSF), PubTator_counts, by.x = "external_gene_name", by.y = "symbol", all.x = T)
# PCSF <- PCSF[order(-PCSF$degree_centrality), ]
# PCSF <- merge(PCSF, druggability, by = "external_gene_name", all.x = T)
# PCSF <- PCSF[complete.cases(PCSF), ]

HHnet_enrich <- merge.data.table(as.data.table(HHnet_enrich), PubTator_counts, by.x = "external_gene_name", by.y = "symbol", all.x = T)
HHnet_enrich <- HHnet_enrich[order(-HHnet_enrich$degree), ]
HHnet_enrich <- merge(HHnet_enrich, druggability, by = "external_gene_name", all.x = T)
temp <- HHnet_enrich[is.na(HHnet_enrich$druggability), ]
HHnet_enrich <- HHnet_enrich[complete.cases(HHnet_enrich), ]

# normalise data
# PCSF$counts_norm <- log10(PCSF$counts)
# PCSF[, c(6:10,17)] <- lapply(PCSF[, c(6:10,17)], function(x) {
#   (x - min(x)) / (max(x) - min(x))
# })

HHnet_enrich$counts_norm <- log10(HHnet_enrich$counts)
HHnet_enrich[, c(5:9,15)] <- lapply(HHnet_enrich[, c(5:9,15)], function(x) {
  (x - min(x)) / (max(x) - min(x))
})


## Average ranking system
avg_rank_sensitivity <- function(input_data,
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


RS_HHnet_enrich <- avg_rank_sensitivity(HHnet_enrich,
                                        features = list(
                                           betweenness = "betweenness",
                                           centrality = "degree_centrality",
                                           druggability = "highest_score",
                                           eigen_centrality = "eigen_centrality",
                                           closeness = "closeness",
                                           page_rank = "page_rank"), 
                                        step = 0.1)

rownames(RS_HHnet_enrich) <- NULL

write.csv(RS_HHnet_enrich, "results/avg_RS_HHnetEnrich.csv", row.names = F)

# run avg rank with PCSF
# RS_PCSF <- avg_rank_sensitivity(PCSF,
#                                 features = list(
#                                   betweenness = "betweenness",
#                                   centrality = "degree_centrality",
#                                   druggability = "highest_score",
#                                   eigen_centrality = "eigen_centrality",
#                                   closeness = "closeness",
#                                   page_rank = "page_rank"), 
#                                 step = 0.1)
# 
# rownames(RS_PCSF) <- NULL
# 
# write.csv(RS_PCSF_new, "results/avg_RS_PCSF.csv", row.names = F)



# analyse overlap with ML results
RS_HHnet_enrich <- read.csv("results/avg_RS_HHnetEnrich.csv")
load("~/OneDrive - RMIT University/PhD/large_git_files/ML/ML_bagging_more_features(100itr).RData")

predicted_targets <- model_prediction_results$feature_data_scores_appended$Protein[model_prediction_results$feature_data_scores_appended$Prediction_Score_rf >= 0.5]
feature_data_scores_appended <- model_prediction_results$feature_data_scores_appended

# annotate with uniprot ids
RS_HHnet_enrich <- id_annot(ensembl, RS_HHnet_enrich, input_type = "external_gene_name", convert_to = "uniprot_gn_id")

RS_HHnet_enrich <- merge(RS_HHnet_enrich, feature_data_scores_appended[, c(1,105:108)], by.x = "uniprot_gn_id", by.y = "Protein", all.x = T)
RS_HHnet_enrich <- merge.data.table(RS_HHnet_enrich, PubTator_counts, by.x = "external_gene_name", by.y = "symbol", all.x = T)
RS_HHnet_enrich <- RS_HHnet_enrich[order(RS_HHnet_enrich$avg_rank), ]
rownames(RS_HHnet_enrich) <- NULL

write.csv(RS_HHnet_enrich, "results/HHnet_RS_ML_overlap.csv", row.names = F)
RS_HHnet_enrich <- read.csv("results/HHnet_RS_ML_overlap.csv")

venn_data <- list(
  HHnetEnrich = RS_HHnet_enrich$uniprot_gn_id,
  #PCSF = PCSF_RS$uniprot_gn_id,
  `ML predicted` = predicted_targets)

library(venn)
library(RColorBrewer)
venn(venn_data, 
     ellipse = T, 
     zcolor = brewer.pal(length(venn_data), name = "Dark2"),
     box = FALSE,
     ilabels = "counts",
     sncs = 2,
     ilcs = 2)



# format a nice table
druggability <- read.csv("../Druggability_analysis/data_general/druggability_scores_annot.csv")
temp <- merge(RS_HHnet_enrich[,c(1:4,6,9)], druggability, by = c("external_gene_name", "uniprot_gn_id"), all.x = TRUE)
temp <- temp[order(temp$avg_rank), ]
rownames(temp) <- NULL
temp$CP_score <- round(temp$CP_score, digits = 3)
temp$`FP/CP` <- paste(temp$druggability, temp$CP_score, sep = "/")
temp <- temp[, -c(8:13)]
temp <- temp[,c("external_gene_name", "uniprot_gn_id", "description", "avg_rank", 
                "rank_variance", "Prediction_Score_rf", "counts", "FP/CP", "source")]
rm(druggability)

# subset genes with less than 5000 citations
temp2 <- temp[temp$counts < 5000, ]

# predicted by RF
temp3 <- temp[temp$Prediction_Score_rf > 0.5, ]

# identify genes with existing/(pre)clinical drugs
load("../Druggability_analysis/data_general/TTD_master.RData")
temp3 <- TTD_master[TTD_master$all_target_genes %in% temp$external_gene_name, ]
temp3 <- temp3[order(match(temp3$all_target_genes, temp$external_gene_name)), ]





