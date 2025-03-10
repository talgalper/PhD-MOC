

target_weight_sensitivity <- function(input_data, target_genes,
                                      features = list(
                                        betweenness = "betweenness",
                                        centrality = "degree_centrality",
                                        druggability = "highest_score",
                                        eigen_centrality = "eigen_centrality",
                                        closeness = "closeness",
                                        citation = "counts_norm"
                                      ),
                                      negative_features = c("citation"),
                                      step = 0.1) {
  library(progress)
  
  # Get feature names and corresponding weight names (e.g. "betweenness" becomes "betweenness_w")
  feature_names <- names(features)
  weight_names <- paste0(feature_names, "_w")
  
  # 1. Generate all weight combinations for selected features
  grid_list <- lapply(feature_names, function(x) seq(0, 1, step))
  names(grid_list) <- weight_names
  all_combos <- expand.grid(grid_list)
  
  # 2. Retain only those combinations that sum to 1 (within numerical tolerance)
  valid_combos <- subset(all_combos, abs(rowSums(all_combos) - 1) < 1e-9)
  
  # Set up a progress bar
  pb <- progress_bar$new(
    format = "[:bar] :current/:total (:percent) eta: :eta",
    total = nrow(valid_combos)
  )
  
  # Preallocate list to store results
  results_list <- vector("list", nrow(valid_combos))
  
  # Iterate over each valid weight combination
  for(i in seq_len(nrow(valid_combos))) {
    pb$tick()
    current_weights <- valid_combos[i, ]
    
    # Compute the combined score dynamically for each gene
    combined_score <- rep(0, nrow(input_data))
    for(feature in feature_names) {
      weight_col <- paste0(feature, "_w")
      sign_factor <- ifelse(feature %in% negative_features, -1, 1)
      combined_score <- combined_score +
        sign_factor * current_weights[[weight_col]] * input_data[[ features[[feature]] ]]
    }
    
    # Compute ranking: higher combined score gets a better (lower) rank
    order_indices <- order(-combined_score)
    ranks <- integer(length(combined_score))
    ranks[order_indices] <- seq_along(combined_score)
    
    # For each target gene, retrieve its rank (if found)
    target_ranks <- sapply(target_genes, function(gene) {
      idx <- which(input_data$external_gene_name == gene)
      if(length(idx) > 0) {
        return(ranks[idx])
      } else {
        return(NA)
      }
    })
    
    # Compute a summary statistic (average rank across target genes)
    avg_rank <- mean(target_ranks, na.rm = TRUE)
    
    # Combine the weight settings and target gene ranks into one result row
    result_row <- c(
      as.list(current_weights),
      setNames(as.list(target_ranks), paste0("rank_", target_genes)),
      avg_rank = avg_rank
    )
    
    results_list[[i]] <- as.data.frame(result_row, stringsAsFactors = FALSE)
  }
  
  # Combine all individual results into a single data frame
  results_df <- do.call(rbind, results_list)
  return(results_df)
}



HHnet <- read.csv("latest_run/OG_rank/HHnet_input.csv", row.names = 1)
HHnet_enrich <- read.csv("~/Desktop/HHnetEnrich_pageRank.csv")
PCSF <- read.csv("latest_run/OG_rank/PCSF_input.csv", row.names = 1)

colnames(HHnet)[5] <- "degree_centrality"
colnames(HHnet_enrich)[5] <- "degree_centrality"

targets <- read.csv("../Druggability_analysis/data_general/target_all_dbs.csv")
targets <- targets[, c(2:4)]
targets <- targets[!duplicated(targets$ensembl_gene_id), ]

# broader set of known breast targets
load("~/Desktop/temp.RData")
TTD_master <- TTD_master[grep("breast", TTD_master$INDICATION, ignore.case = T), ]
TTD_master <- TTD_master[TTD_master$HighestClinicalStatus == "Approved", ]
TTD_master <- unique(TTD_master$all_target_genes)
TTD_master <- na.omit(TTD_master)
TTD_master <- strsplit(TTD_master, ";")
TTD_master <- unlist(TTD_master)
TTD_master <- strsplit(TTD_master, "/")
TTD_master <- unlist(TTD_master)
TTD_master <- unique(TTD_master)
TTD_master <- trimws(TTD_master)
table(TTD_master %in% HHnet_enrich$external_gene_name)
table(targets$drugBank_target %in% TTD_master)



# normalise data
PCSF$counts_norm <- log10(PCSF$counts)
PCSF[, c(6:9,16)] <- lapply(PCSF[, c(6:9,16)], function(x) {
  (x - min(x)) / (max(x) - min(x))
})

HHnet$counts_norm <- log10(HHnet$counts)
HHnet[, c(5:8,15)] <- lapply(HHnet[, c(5:8,15)], function(x) {
  (x - min(x)) / (max(x) - min(x))
})

HHnet_enrich$counts_norm <- log10(HHnet_enrich$counts)
HHnet_enrich[, c(5:8,15)] <- lapply(HHnet_enrich[, c(5:8,15)], function(x) {
  (x - min(x)) / (max(x) - min(x))
})


TWS_HHnet_enrich <- target_weight_sensitivity(HHnet_enrich, unique(targets$drugBank_target))
TWS_PCSF <- target_weight_sensitivity(PCSF, unique(targets$drugBank_target))

# run again without citation scores
TWS_HHnet_enrich_noCite <- target_weight_sensitivity(HHnet_enrich, unique(targets$drugBank_target),
                                              features = list(
                                                betweenness = "betweenness",
                                                centrality = "degree_centrality",
                                                druggability = "highest_score",
                                                eigen_centrality = "eigen_centrality",
                                                closeness = "closeness"), 
                                              step = 0.1)

TWS_PCSF_noCite <- target_weight_sensitivity(PCSF, unique(targets$drugBank_target),
                                      features = list(
                                        betweenness = "betweenness",
                                        centrality = "degree_centrality",
                                        druggability = "highest_score",
                                        eigen_centrality = "eigen_centrality",
                                        closeness = "closeness"), 
                                      step = 0.1)


# run again with page_rank
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

# Assuming 'results_df' is the output from target_weight_sensitivity()
features <- c("betweenness_w", "centrality_w", "druggability_w", "eigen_centrality_w", "closeness_w", "page_rank_w")


# Initialize an empty list to store the results
results_list <- lapply(features, function(feature) {
  # Run the Spearman correlation test
  test <- cor.test(TWS_HHnet_enrich_pageRank[[feature]], TWS_HHnet_enrich_pageRank$avg_rank, method = "spearman")
  
  # Extract the correlation estimate, p-value, and confidence interval (if available)
  conf_int <- if(!is.null(test$conf.int)) test$conf.int else c(NA, NA)
  
  data.frame(
    feature = feature,
    spearman_rho = test$estimate,
    p_value = test$p.value,
    conf_low = conf_int[1],
    conf_high = conf_int[2],
    stringsAsFactors = FALSE
  )
})

# Combine the results into a single data frame
HHnet_results_df <- do.call(rbind, results_list)


# Repeat for PCSF
results_list <- lapply(features, function(feature) {
  # Run the Spearman correlation test
  test <- cor.test(TWS_PCSF_pageRank[[feature]], TWS_PCSF_pageRank$avg_rank, method = "spearman")
  
  # Extract the correlation estimate, p-value, and confidence interval (if available)
  conf_int <- if(!is.null(test$conf.int)) test$conf.int else c(NA, NA)
  
  data.frame(
    feature = feature,
    spearman_rho = test$estimate,
    p_value = test$p.value,
    conf_low = conf_int[1],
    conf_high = conf_int[2],
    stringsAsFactors = FALSE
  )
})

# Combine the results into a single data frame
PCSF_results_df <- do.call(rbind, results_list)








top_results <- TWS_HHnet_enrich_noCite[order(TWS_HHnet_enrich_noCite$avg_rank), ]
top_results <- head(top_results, n = 10)
top_results$combo <- paste0("C_", top_results$centrality_w, 
                           "/eC_", top_results$eigen_centrality_w, 
                           "/B_", top_results$betweeness_w,
                           "/Cl_", top_results$closeness_w,
                           "/D_", top_results$druggability_w,
                           "/ci_", top_results$citation_w)
ggplot(top_results, aes(x = reorder(combo, avg_rank), y = avg_rank)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top 10 Weight Combinations",
       x = "Weight Combination",
       y = "Average Rank")

