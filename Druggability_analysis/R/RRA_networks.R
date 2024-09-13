library(igraph)
library(tidyverse)
library(PCSF)

druggability <- read.csv("data_general/druggability_source.csv", row.names = 1)

targets <- read.csv("../Druggability_analysis/data_general/target_all_dbs.csv")
targets <- unique(targets$drugBank_target)

HHnet <- read_graph("../Hierarchical_HotNet/BRCA/results/hhnet_cluster1_netNeighs.graphml", format = "graphml")
HHnet <- as.undirected(HHnet)

druggability_scores <- druggability$highest_score[match(V(HHnet)$name, druggability$external_gene_name)]
druggability_scores[is.na(druggability_scores)] <- 0

HHnet <- set_vertex_attr(HHnet, name = "druggability", value = druggability_scores)

#V(HHnet)$size <- scales::rescale(V(HHnet)$druggability, c(1, 5)) # resize based on druggability, this is cooked
V(HHnet)$color <- ifelse(V(HHnet)$name %in% targets, "blue", V(HHnet)$color) # colour known targets
vertex_labels <- ifelse(V(HHnet)$name %in% targets, V(HHnet)$name, NA) # plot only known target labels

plot(HHnet, asp = 0, vertex.label = vertex_labels, edge.arrow.size = 0.3, vertex.label.dist = 1) 


df_subnetNeighs <- data.frame(display.name = V(HHnet)$`display name`,
                              degree = degree(HHnet),
                              betweenness = betweenness(HHnet),
                              closeness = closeness(HHnet),
                              eigen_centrality = eigen_centrality(HHnet)$vector,
                              drug_poc_score = V(HHnet)$druggability,
                              source = V(HHnet)$color)

df_subnetNeighs$source <- ifelse(df_subnetNeighs$source == "tomato", "subnet", "STRING")
df_subnetNeighs$source <- ifelse(df_subnetNeighs$source == "blue", "known_target", df_subnetNeighs$source)
df_subnetNeighs <- df_subnetNeighs[order(-df_subnetNeighs$degree), ]
rownames(df_subnetNeighs) <- NULL

rm(HHnet)
gc()

load("../BRCA_pipe/latest_run/RData/PCSF_subnet_FULL.RData")


# extract cluster data
PCSF_df <- data.frame(gene_id = V(subnet)$name,
                      degree = degree(subnet),
                      betweenness = betweenness(subnet),
                      closeness = closeness(subnet),
                      eigen_centrality = eigen_centrality(subnet)$vector,
                      prize = V(subnet)$prize,
                      source = V(subnet)$type)

PCSF_df <- PCSF_df[order(-PCSF_df$degree), ]
rownames(PCSF_df) <- NULL

rm(subnet)
gc()


df <- df_subnetNeighs[df_subnetNeighs$display.name %in% targets, ]
df <- PCSF_df[PCSF_df$gene_id %in% targets, ]


library(RobustRankAggreg)
RRA_HHnet <- aggregateRanks(list(degree = df_subnetNeighs$display.name[order(-df_subnetNeighs$degree)],
                           betweenness = df_subnetNeighs$display.name[order(-df_subnetNeighs$betweenness)],
                           closeness = df_subnetNeighs$display.name[order(-df_subnetNeighs$closeness)],
                           eigen_centrality = df_subnetNeighs$display.name[order(-df_subnetNeighs$eigen_centrality)]))
rownames(RRA_HHnet) <- NULL
RRA_HHnet <- RRA_HHnet[RRA_HHnet$Name %in% targets, ]


RRA_PCSF <- aggregateRanks(list(degree = PCSF_df$gene_id[order(-PCSF_df$degree)],
                           betweenness = PCSF_df$gene_id[order(-PCSF_df$betweenness)],
                           closeness = PCSF_df$gene_id[order(-PCSF_df$closeness)],
                           eigen_centrality = PCSF_df$gene_id[order(-PCSF_df$eigen_centrality)]))
rownames(RRA_PCSF) <- NULL
RRA_PCSF <- RRA_PCSF[RRA_PCSF$Name %in% targets, ]





### test all combinations of features in RRA ###
library(RobustRankAggreg)

# Define all metrics
metrics <- list(
  degree = df_subnetNeighs$display.name[order(-df_subnetNeighs$degree)],
  betweenness = df_subnetNeighs$display.name[order(-df_subnetNeighs$betweenness)],
  closeness = df_subnetNeighs$display.name[order(-df_subnetNeighs$closeness)],
  eigen_centrality = df_subnetNeighs$display.name[order(-df_subnetNeighs$eigen_centrality)]
)

metrics <- list(
  degree = PCSF_df$gene_id[order(-PCSF_df$degree)],
  betweenness = PCSF_df$gene_id[order(-PCSF_df$betweenness)],
  closeness = PCSF_df$gene_id[order(-PCSF_df$closeness)],
  eigen_centrality = PCSF_df$gene_id[order(-PCSF_df$eigen_centrality)]
)

# with druggability scores
temp <- druggability[, c(2,8)]
PCSF_df2 <- merge(PCSF_df, temp, by.x = "gene_id", by.y = "external_gene_name")

metrics <- list(
  degree = df_subnetNeighs2$display.name[order(-df_subnetNeighs2$degree)],
  betweenness = df_subnetNeighs2$display.name[order(-df_subnetNeighs2$betweenness)],
  closeness = df_subnetNeighs2$display.name[order(-df_subnetNeighs2$closeness)],
  eigen_centrality = df_subnetNeighs2$display.name[order(-df_subnetNeighs2$eigen_centrality)],
  druggability = df_subnetNeighs2$display.name[order(-df_subnetNeighs2$drug_poc_score)]
)

metrics <- list(
  degree = PCSF_df2$gene_id[order(-PCSF_df2$degree)],
  betweenness = PCSF_df2$gene_id[order(-PCSF_df2$betweenness)],
  closeness = PCSF_df2$gene_id[order(-PCSF_df2$closeness)],
  eigen_centrality = PCSF_df2$gene_id[order(-PCSF_df2$eigen_centrality)],
  druggability = PCSF_df2$gene_id[order(-PCSF_df2$highest_score)]
)



# Function to compute RRA for a given set of metrics and count top 100 targets
evaluate_combination <- function(selected_metrics, targets) {
  RRA <- aggregateRanks(selected_metrics)
  rownames(RRA) <- NULL
  RRA <- RRA[RRA$Name %in% targets, ]
  RRA$rank <- as.integer(rownames(RRA))
  
  # Count how many targets are in the top 100
  sum(RRA$rank <= 100)
}

# Get all combinations of metrics
metric_combinations <- unlist(lapply(1:length(metrics), function(x) combn(names(metrics), x, simplify = FALSE)), recursive = FALSE)

# Evaluate each combination and store results
results <- sapply(metric_combinations, function(metric_names) {
  selected_metrics <- metrics[metric_names]
  evaluate_combination(selected_metrics, targets)
})

# Find the best combination
max_score <- max(results)
best_combinations_indices <- which(results == max_score)
best_combinations <- metric_combinations[best_combinations_indices]
best_combinations