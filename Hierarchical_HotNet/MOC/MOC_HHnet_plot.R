library(tidyverse)
library(igraph)

STRING_net <- read_graph("../../../../OneDrive - RMIT University/PhD/large_git_files/HHnet/STRING_hsa_physical_network.graphml", format = "graphml") # ubuntu
STRING_net <- as.undirected(STRING_net)

hh_results <- read_lines("MOC/results/clusters_STRING_MOC_logFC_scores_abs.tsv", skip = 7)
hh_results <- str_split(hh_results, pattern = "\t")

#sum(hh_results[[1]] %in% V(STRING_net)$`display name`) 

# colour cluster members
V(STRING_net)$color <- ifelse(V(STRING_net)$`display name` %in% hh_results[[1]], "tomato", "white")

#idx <- unique(c(
#  which(V(STRING_net)$color == "tomato"),
#  seq(vcount(STRING_net))
#))

# too much going on when trying to plot everything
#STRING_net_plot <- permute(STRING_net, rev(idx))
#STRING_net_plot <- simplify(STRING_net_plot)
#STRING_net_plot <- induced.subgraph(STRING_net_plot, V(STRING_net_plot)[degree(STRING_net_plot) > 0]) # remove nodes with degree = 0
#plot(STRING_net_plot, asp = 0, vertex.label = NA, vertex.size = 0.5, edge.arrow.size = 0.1)



# subnetwork cluster
subnet <- induced_subgraph(graph = STRING_net, V(STRING_net)$`display name` %in% hh_results[[1]])
plot.igraph(subnet, asp = 0, vertex.size = 2, edge.arrow.size = 0.3, vertex.label.dist = 1, vertex.label = V(subnet)$`display name`)

# add direct neighbors of cluster nodes
neighbs <- neighborhood(
  STRING_net,
  order = 1,
  nodes = hh_results[[1]]
)
neighbs <- flatten(neighbs)
neighbs <- names(neighbs)
neighbs <- unique(neighbs)

clust1_net <- induced_subgraph(STRING_net, V(STRING_net)[V(STRING_net)$`display name` %in% neighbs])
V(clust1_net)$size <- scales::rescale(degree(clust1_net), c(1, 7))
#clust1_net <- simplify(clust1_net)
plot.igraph(clust1_net, asp = 0, vertex.label = NA, edge.arrow.size = 0.3)

# subset edges with experimental evidence
#clust1_net_filt <- delete_edges(clust1_net, E(clust1_net)[!is.na(E(clust1_net)$`stringdb::experiments`) & E(clust1_net)$`stringdb::experiments` >= 0.4])
#plot(clust1_net_filt, asp = 0, vertex.label = NA, vertex.size = 2, edge.arrow.size = 0.3)

# delete the ENSP IDs for later, bites you in the ass with MCODE
#subnet <- delete_vertex_attr(subnet, "shared name")
#clust1_net <- delete_vertex_attr(clust1_net, "shared name")

write_graph(clust1_net, "MOC/results/hhnet_cluster1_netNeighs.graphml", format = "graphml")
write_graph(subnet, "MOC/results/hhnet_cluster1_net.graphml", format = "graphml")


df <- data.frame(degree = degree(subnet),
                 betweenness = betweenness(subnet),
                 source = V(subnet)$color)
df$source <- ifelse(df$source == "tomato", "subnet", "STRING")
df <- rownames_to_column(df)
df <- df[order(-df$degree), ]


df <- data.frame(degree = degree(clust1_net),
                 betweenness = betweenness(clust1_net),
                 source = V(clust1_net)$color)
df$source <- ifelse(df$source == "tomato", "subnet", "STRING")
df <- rownames_to_column(df)
df <- df[order(-df$degree), ]







# analyse MCODE clusters
MCODE_cluster1 <- read_lines("results/MCODE_cluster1.txt", skip = 9)
MCODE_cluster1 <- str_split(MCODE_cluster1, "\t")
MCODE_cluster1 <- do.call(rbind, MCODE_cluster1)
MCODE_cluster1 <- as.data.frame(MCODE_cluster1)
colnames(MCODE_cluster1) <- MCODE_cluster1[1, ]
MCODE_cluster1 <- MCODE_cluster1[-1, ]
MCODE_clusters <- str_split(MCODE_cluster1$`Node IDs`, ", ")
names(MCODE_clusters) <- 1:length(MCODE_clusters)

# analyse MCODE clusters with neighbors
MCODE_cluster1_neighs <- read_lines("results/MCODE_cluster1_neighs.txt", skip = 9)
MCODE_cluster1_neighs <- str_split(MCODE_cluster1_neighs, "\t")
MCODE_cluster1_neighs <- do.call(rbind, MCODE_cluster1_neighs)
MCODE_cluster1_neighs <- as.data.frame(MCODE_cluster1_neighs)
colnames(MCODE_cluster1_neighs) <- MCODE_cluster1_neighs[1, ]
MCODE_cluster1_neighs <- MCODE_cluster1_neighs[-1, ]
MCODE_neigh_clusters <- str_split(MCODE_cluster1_neighs$`Node IDs`, ", ")
names(MCODE_neigh_clusters) <- 1:length(MCODE_neigh_clusters)

# Loop through each cluster and assign the corresponding color
cluster_colors <- grDevices::rainbow(length(MCODE_clusters))
names(cluster_colors) <- names(MCODE_clusters)
V(subnet)$color <- "grey"  # Default color for nodes not in any MCODE cluster

for (i in seq_along(MCODE_clusters)) {
  V(subnet)$color[V(subnet)$name %in% MCODE_clusters[[i]]] <- cluster_colors[i]
}
plot(subnet, asp = 0, vertex.label = NA, vertex.size = 2, edge.arrow.size = 0.3)


library(clusterProfiler)
library(org.Hs.eg.db)
library(progress)

# Go on whole networks
GO <- enrichGO(V(subnet)$name, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP")
GO_formatted <- GO@result
GO2 <- enrichGO(V(clust1_net)$name, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP")
GO2_formatted <- GO2@result


pb <- progress_bar$new(total = length(MCODE_clusters))

# run GO enrichment
# turn this into a function at some point
GO_clust1 <- list()
for (cluster in MCODE_clusters) {
  GO <- enrichGO(cluster, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP")
  
  GO_clust1 <- append(GO_clust1, GO)
  
  rm(cluster, GO)
  pb$tick()
}
names(GO_clust1) <- 1:length(GO_clust1)


GO_clust1_formatted <- data.frame()
for (i in seq_along(GO_clust1)) {
  cluster <- GO_clust1[[i]]
  cluster_name <- names(GO_clust1)[i]
  
  result <- cluster@result
  result_top <- head(result, 5)
  result_top$cluster <- rep(cluster_name, nrow(result_top))
  result_top$`-log(p.adjust)` <- -log(result_top$p.adjust)
  
  GO_clust1_formatted <- rbind(GO_clust1_formatted, result_top)
  
  rm(cluster, i, result, result_top, cluster_name)
}

# Function to convert GeneRatio to numeric
convert_gene_ratio <- function(gene_ratio) {
  # Split the string by "/"
  parts <- strsplit(gene_ratio, "/")[[1]]
  # Convert the parts to numeric and calculate the ratio
  ratio <- as.numeric(parts[1]) / as.numeric(parts[2])
  return(ratio)
}
GO_clust1_formatted$GeneRatio.num <- sapply(GO_clust1_formatted$GeneRatio, convert_gene_ratio)

save(GO_clust1_formatted, file = "results/RData/GO_clust1_formatted.RData")

# too messy
#ggplot(data = GO_clust1_formatted, aes(x = cluster, y = Description, 
#                                color = `p.adjust`, size = GeneRatio.num)) + 
#  geom_point() +
#  scale_color_gradient(low = "red", high = "blue") +
#  theme_bw() + 
#  ylab("") + 
#  xlab("") + 
#  labs(size = "Gene Ratio") +
#  ggtitle("GO enrichment analysis (BP)")

