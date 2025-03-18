library(tidyverse)
library(progress)
library(biomaRt)

rank_sensitivity <- function(input_data, ensembl, topx,
                             # 'features' is a named list where the name is used to build weight columns
                             # and the value is the corresponding column name in input_data.
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
  
  # Extract feature names and generate corresponding weight column names
  feature_names <- names(features)
  weight_names <- paste0(feature_names, "_w")
  
  # Dynamically generate all weight combinations for selected features
  grid_list <- lapply(seq_along(feature_names), function(i) seq(0, 1, step))
  names(grid_list) <- weight_names
  all_combos <- expand.grid(grid_list)
  
  # Keep only the combinations where the weights sum to 1 (within tolerance)
  valid_combos <- subset(all_combos, abs(rowSums(all_combos) - 1) < 1e-9)
  
  # Create a progress bar over valid_combos only
  pb <- progress_bar$new(
    format = "[:bar] :current/:total (:percent) eta: :eta", 
    total = nrow(valid_combos)
  )
  
  # Preallocate list for results
  results_list <- vector("list", length = nrow(valid_combos))
  
  # Loop over each valid combination
  for(i in seq_len(nrow(valid_combos))) {
    pb$tick()
    current_weights <- valid_combos[i, ]
    
    # Compute combined score dynamically based on selected features.
    # For each feature, multiply the weight by the corresponding input_data column.
    # Use a negative sign for features specified in negative_features.
    combined_score <- rep(0, nrow(input_data))
    for(feature in feature_names) {
      weight_col <- paste0(feature, "_w")
      sign_factor <- ifelse(feature %in% negative_features, -1, 1)
      combined_score <- combined_score +
        sign_factor * current_weights[[weight_col]] * input_data[[ features[[feature]] ]]
    }
    
    # Extract the top 'topx' indices (using partial sorting)
    top_idx <- order(combined_score, decreasing = TRUE)[1:topx]
    top_genes <- input_data$external_gene_name[top_idx]
    
    # Build a result row: include the weight settings and the top genes as a comma‐separated string.
    results_list[[i]] <- data.frame(
      # Dynamically add the weight values to the result row
      setNames(as.list(current_weights), weight_names),
      top_genes = paste(top_genes, collapse = ', '),
      stringsAsFactors = FALSE
    )
  }
  
  # Combine results into one data frame
  sensitivity_results <- do.call(rbind, results_list)
  
  # Split the top genes into lists for further counting
  sensitivity_results$top_genes_list <- strsplit(sensitivity_results$top_genes, ', ')
  
  # Count gene occurrences across all combinations
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
  
  # For each weight, count the appearances of each gene when that weight is 0
  for (w in weight_names) {
    other_weights <- setdiff(weight_names, w)
    subset_df <- sensitivity_results[
      sensitivity_results[[w]] == 0 & rowSums(sensitivity_results[, other_weights] == 0) == 0, ]
    
    counts_col <- sapply(gene_counts$gene, function(gene) {
      sum(sapply(subset_df$top_genes_list, function(lst) gene %in% lst))
    })
    
    new_col_name <- paste0("counts_when_", sub("_w$", "", w), "_0")
    gene_counts[[new_col_name]] <- counts_col
  }
  rownames(gene_counts) <- NULL
  
  # Retrieve gene descriptions from Ensembl
  gene_description <- getBM(
    attributes = c("external_gene_name", "description"), 
    filters = "external_gene_name", 
    values = gene_counts$gene,
    mart = ensembl
  )
  gene_description$description <- gsub("\\s*\\[.*?\\]", "", gene_description$description)
  
  final_gene_counts <- merge(gene_description, gene_counts, 
                             by.x = "external_gene_name", by.y = "gene", all.y = TRUE)
  final_gene_counts <- final_gene_counts[order(-final_gene_counts$count), ]
  rownames(final_gene_counts) <- NULL
  
  return(final_gene_counts)
}




ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

HHnet <- read.csv("latest_run/OG_rank/HHnet_input.csv", row.names = 1)
HHnet_enrich <- read.csv("latest_run/OG_rank/HHnet_enrich_input.csv", row.names = 1)
PCSF <- read.csv("latest_run/OG_rank/PCSF_input.csv", row.names = 1)

# column names must match
colnames(HHnet)[5] <- "degree_centrality"
colnames(HHnet_enrich)[5] <- "degree_centrality"

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

RS_PCSF <- rank_sensitivity(PCSF, ensembl, topx = 10)
RS_HHnet <- rank_sensitivity(HHnet, ensembl, topx = 10)
RS_HHnet_enrich <- rank_sensitivity(HHnet_enrich, ensembl, topx = 10)

write.csv(RS_PCSF, "latest_run/OG_rank/PCSF_rank_sensitivity_top10.csv", row.names = F)
write.csv(RS_HHnet, "latest_run/OG_rank/HHnet_rank_sensitivity_top10.csv", row.names = F)
write.csv(RS_HHnet_enrich, "latest_run/OG_rank/HHnetEnrich_rank_sensitivity_top10.csv", row.names = F)

RS_PCSF <- read.csv("latest_run/OG_rank/PCSF_rank_sensitivity_top10.csv")
RS_HHnet <- read.csv("latest_run/OG_rank/HHnet_rank_sensitivity_top10.csv")
RS_HHnet_enrich <- read.csv("latest_run/OG_rank/HHnetEnrich_rank_sensitivity_top10.csv")

# number of targets in RS results
targets <- read.csv("../Druggability_analysis/data_general/target_all_dbs.csv")
targets <- targets[, c(2,4)]
targets <- targets[!duplicated(targets$ensembl_gene_id), ]

temp <- RS_HHnet_enrich[RS_HHnet_enrich$external_gene_name %in% targets$drugBank_target, ]
temp2 <- RS_PCSF[RS_PCSF$external_gene_name %in% targets$drugBank_target, ]

# intersect between RS results
library(ggVennDiagram)
p <- ggVennDiagram(list(HHnet = RS_HHnet$external_gene_name,
                        `HHnet Neighs` = RS_HHnet_enrich$external_gene_name,
                        PCSF = RS_PCSF$external_gene_name),
                   label = "count",
                   set_size = 8, label_size = 6)

p + scale_fill_distiller(palette = "Oranges", direction = 1) + coord_equal(clip = "off") +
  theme(
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 20)  
  )


intersect(RS_HHnet$external_gene_name, intersect(RS_HHnet_enrich$external_gene_name, RS_PCSF$external_gene_name))
PubTator_counts[symbol %in% c("FHIP1B", "GPR156", "ACAN", "HAO2")]


i <- intersect(RS_HHnet_enrich$external_gene_name, RS_PCSF$external_gene_name)
PubTator_counts[symbol %in% i]$symbol

i <- intersect(RS_HHnet_enrich_no_cit$external_gene_name, RS_PCSF_no_cit$external_gene_name)
PubTator_counts[symbol %in% i]



# attempt with no citation score
RS_PCSF_noCite <- rank_sensitivity(PCSF, ensembl, topx = 10, 
                                   features = list(
                                     betweenness = "betweenness",
                                     centrality = "degree_centrality",
                                     druggability = "highest_score",
                                     eigen_centrality = "eigen_centrality",
                                     closeness = "closeness"), 
                                   step = 0.1)

RS_HHnet_noCite <- rank_sensitivity(HHnet, ensembl, topx = 10, 
                                    features = list(
                                      betweenness = "betweenness",
                                      centrality = "degree_centrality",
                                      druggability = "highest_score",
                                      eigen_centrality = "eigen_centrality",
                                      closeness = "closeness"), 
                                    step = 0.1)

RS_HHnet_enrich_noCite <- rank_sensitivity(HHnet_enrich, ensembl, topx = 10, 
                                           features = list(
                                             betweenness = "betweenness",
                                             centrality = "degree_centrality",
                                             druggability = "highest_score",
                                             eigen_centrality = "eigen_centrality",
                                             closeness = "closeness"), 
                                           step = 0.1)

write.csv(RS_HHnet_enrich_noCite, "latest_run/OG_rank/HHnet_rank_sensitivity_noCite.csv", row.names = F)
write.csv(RS_PCSF_noCite, "latest_run/OG_rank/PCSF_rank_sensitivity_noCite.csv", row.names = F)

temp <- RS_HHnet_enrich_noCite[RS_HHnet_enrich_noCite$external_gene_name %in% targets$drugBank_target, ]
temp2 <- RS_PCSF_noCite[RS_PCSF_noCite$external_gene_name %in% targets$drugBank_target, ]


# attempt with increased RS topx
RS_PCSF_noCite_60 <- rank_sensitivity(PCSF, ensembl, topx = 60, 
                                   features = list(
                                     betweenness = "betweenness",
                                     centrality = "degree_centrality",
                                     druggability = "highest_score",
                                     eigen_centrality = "eigen_centrality",
                                     closeness = "closeness"), 
                                   step = 0.1)

RS_HHnet_noCite_60 <- rank_sensitivity(HHnet, ensembl, topx = 60, 
                                    features = list(
                                      betweenness = "betweenness",
                                      centrality = "degree_centrality",
                                      druggability = "highest_score",
                                      eigen_centrality = "eigen_centrality",
                                      closeness = "closeness"), 
                                    step = 0.1)

RS_HHnet_enrich_noCite <- rank_sensitivity(HHnet_enrich, ensembl, topx = 10, 
                                           features = list(
                                             betweenness = "betweenness",
                                             centrality = "degree_centrality",
                                             druggability = "highest_score",
                                             eigen_centrality = "eigen_centrality",
                                             closeness = "closeness"), 
                                           step = 0.1)

temp <- RS_HHnet_enrich_noCite_60[RS_HHnet_enrich_noCite_60$external_gene_name %in% targets$drugBank_target, ]
temp2 <- RS_PCSF_noCite_60[RS_PCSF_noCite_60$external_gene_name %in% targets$drugBank_target, ]
