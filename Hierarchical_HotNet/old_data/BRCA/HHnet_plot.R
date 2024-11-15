library(tidyverse)
library(igraph)

STRING_net <- read_graph("~/OneDrive - RMIT University/PhD/large_git_files/HHnet/old/STRING_hsa_physical_network.graphml", format = "graphml")

hh_results <- read_lines("BRCA/results/clusters_STRING_logFC_scores_abs.tsv", skip = 7)
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
plot(subnet, asp = 0, vertex.label = NA, vertex.size = 2, edge.arrow.size = 0.3)

# add direct neighbors of cluster nodes
V(STRING_net)$name <- V(STRING_net)$`display name`
neighbs <- neighborhood(
  STRING_net,
  order = 1,
  nodes = hh_results[[1]]
)
neighbs <- flatten(neighbs)
neighbs <- names(neighbs)
neighbs <- unique(neighbs)

clust1_net <- induced_subgraph(STRING_net, V(STRING_net)[V(STRING_net)$name %in% neighbs])
V(clust1_net)$size <- scales::rescale(degree(clust1_net), c(1, 7))
#clust1_net <- simplify(clust1_net) # pretty sure gets rid of all extra data from STRING/Cytoscape
plot(clust1_net, asp = 0, vertex.label = NA, edge.arrow.size = 0.3)

# subset edges with experimental evidence
#clust1_net_filt <- delete_edges(clust1_net, E(clust1_net)[!is.na(E(clust1_net)$`stringdb::experiments`) & E(clust1_net)$`stringdb::experiments` >= 0.4])
#plot(clust1_net_filt, asp = 0, vertex.label = NA, vertex.size = 2, edge.arrow.size = 0.3)


write_graph(clust1_net, "BRCA/results/hhnet_cluster1_netNeighs.graphml", format = "graphml")
write_graph(subnet, "BRCA/results/hhnet_cluster1_net.graphml", format = "graphml")

clust1_net <- read_graph("BRCA/results/hhnet_cluster1_netNeighs.graphml", format = "graphml")
subnet <- read_graph("BRCA/results/hhnet_cluster1_net.graphml", format = "graphml")


df_subnet <- data.frame(display.name = V(subnet)$`display name`,
                        degree = degree(subnet),
                        betweenness = betweenness(subnet),
                        closeness = closeness(subnet),
                        eigen_centrality = eigen_centrality(subnet)$vector,
                        source = V(subnet)$color)
df_subnet$source <- ifelse(df_subnet$source == "tomato", "subnet", "STRING")
df_subnet <- df_subnet[order(-df_subnet$degree), ]


df_subnetNeighs <- data.frame(display.name = V(clust1_net)$`display name`,
                              degree = degree(clust1_net),
                              betweenness = betweenness(clust1_net),
                              closeness = closeness(clust1_net),
                              eigen_centrality = eigen_centrality(clust1_net)$vector,
                              source = V(clust1_net)$color)
df_subnetNeighs$source <- ifelse(df_subnetNeighs$source == "tomato", "subnet", "STRING")
df_subnetNeighs <- df_subnetNeighs[order(-df_subnetNeighs$degree), ]
rownames(df_subnetNeighs) <- NULL
write.csv(df_subnet, "BRCA/results/df_subnet.csv")
write.csv(df_subnetNeighs, "BRCA/results/df_subnetNeighs.csv")



# network metrics of known BRCA targets
df_subnet <- read.csv("BRCA/results/df_subnet.csv")
df_subnetNeighs <- read.csv("BRCA/results/df_subnetNeighs.csv", row.names = 1)

targets <- read.csv("../Druggability_analysis/data_general/target_all_dbs.csv")
targets <- unique(targets$drugBank_target)

df_subnet <- df_subnet[df_subnet$display.name %in% targets, ]
df_subnetNeighs <- df_subnetNeighs[df_subnetNeighs$display.name %in% targets, ]


# plot the neighbours of the known targets
neighbs <- neighborhood(
  clust1_net,
  order = 1,
  nodes = df_subnetNeighs$display.name
)
neighbs <- flatten(neighbs)
neighbs <- names(neighbs)
neighbs <- unique(neighbs)

temp <- induced_subgraph(clust1_net, V(clust1_net)[V(clust1_net)$name %in% neighbs])
temp <- as.undirected(temp)

V(temp)$color <- ifelse(V(temp)$name %in% df_subnetNeighs$display.name, "blue", V(temp)$color)

plot(temp, asp = 0, vertex.label = NA, edge.arrow.size = 0.3)

temp_experimental <- delete_edges(temp, E(temp)[!is.na(E(temp)$`stringdb::experiments`) & E(temp)$`stringdb::experiments` >= 0.4])
temp_experimental <- induced_subgraph(temp_experimental, V(temp_experimental)[degree(temp_experimental) > 0]) # remove nodes with degree = 0
temp_experimental <- as.undirected(temp_experimental)

plot(temp_experimental, asp = 0, vertex.label = NA, edge.arrow.size = 0.3)


df_temp <- data.frame(display.name = V(temp)$`display name`,
                              degree = degree(temp),
                              betweenness = betweenness(temp),
                              closeness = closeness(temp),
                              eigen_centrality = eigen_centrality(temp)$vector,
                              source = V(temp)$color)

df_tempExperimental <- data.frame(display.name = V(temp_experimental)$`display name`,
                      degree = degree(temp_experimental),
                      betweenness = betweenness(temp_experimental),
                      closeness = closeness(temp_experimental),
                      eigen_centrality = eigen_centrality(temp_experimental)$vector,
                      source = V(temp_experimental)$color)





important_nodes <- which(V(temp)$color == c("tomato", "blue"))
other_nodes <- which(V(temp)$color != c("tomato", "blue"))

# Reorder nodes so that important nodes are plotted last
new_order <- c(other_nodes, important_nodes)

# Plot the graph, with reordered vertices
plot(temp, vertex.label = NA, asp = 0, edge.arrow.size = 0.3,
     vertex.color = V(temp)$color[new_order], layout = layout_with_fr(temp)[new_order, ])



GO <- enrichGO(V(subnet)$name, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP")
GO_formatted <- GO@result
GO2 <- enrichGO(V(clust1_net)$name, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP")
GO2_formatted <- GO2@result





library(RobustRankAggreg)
RRA <- aggregateRanks(list(degree = df_subnetNeighs$display.name[order(-df_subnetNeighs$degree)],
                           betweenness = df_subnetNeighs$display.name[order(-df_subnetNeighs$betweenness)],
                           closeness = df_subnetNeighs$display.name[order(-df_subnetNeighs$closeness)],
                           eigen_centrality = df_subnetNeighs$display.name[order(-df_subnetNeighs$eigen_centrality)]))
rownames(RRA) <- NULL
RRA <- RRA[RRA$Name %in% targets, ]


df_subnetNeighs_edge <- as_data_frame(clust1_net, what = "edges")
df_subnetNeighs_edge$edge_betweenness <- edge_betweenness(clust1_net)

edge_connectivity_values <- sapply(E(clust1_net), function(e) {
  vertices <- ends(clust1_net, e)
  edge_connectivity(clust1_net, source = vertices[1], target = vertices[2])
})

df_subnetNeighs_edge$edge_connectivity <- edge_connectivity_values


# add druggability data
drug_scores <- read.csv("../Druggability_analysis/data_general/druggability_scores_annot.csv")























