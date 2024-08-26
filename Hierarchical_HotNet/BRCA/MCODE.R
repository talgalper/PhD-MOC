


# delete the ENSP IDs for later, bites you in the ass with MCODE
subnet <- delete_vertex_attr(subnet, "shared name")
clust1_net <- delete_vertex_attr(clust1_net, "shared name")

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
GO <- GO@result
GO2 <- enrichGO(V(clust1_net)$name, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP")
GO2 <- GO2@result


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