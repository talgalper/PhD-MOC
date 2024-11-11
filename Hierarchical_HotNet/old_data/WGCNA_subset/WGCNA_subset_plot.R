library(tidyverse)
library(igraph)

STRING_net <- read_graph("../../../../Desktop/hierarchical-hotnet/STRING_hsa_physical_network.graphml", format = "graphml") # ubuntu
STRING_net <- read_graph("../../../../OneDrive - RMIT University/PhD/large_git_files/HHnet/STRING_hsa_physical_network.graphml", format = "graphml") # MAC
STRING_net <- as.undirected(STRING_net)

hh_results <- read_lines("WGCNA_subset/results/clusters_STRING_MOC_logFC_scores_abs.tsv", skip = 7)
hh_results <- str_split(hh_results, pattern = "\t")

#sum(hh_results[[1]] %in% V(STRING_net)$`display name`) 

# colour cluster members
V(STRING_net)$color <- ifelse(V(STRING_net)$`display name` %in% hh_results[[1]], "tomato", "white")

idx <- unique(c(
  which(V(STRING_net)$color == "tomato"),
  seq(vcount(STRING_net))
))

# too much going on when trying to plot everything
#STRING_net_plot <- permute(STRING_net, rev(idx))
#STRING_net_plot <- simplify(STRING_net_plot)
#STRING_net_plot <- induced.subgraph(STRING_net_plot, V(STRING_net_plot)[degree(STRING_net_plot) > 0]) # remove nodes with degree = 0
#plot(STRING_net_plot, asp = 0, vertex.label = NA, vertex.size = 0.5, edge.arrow.size = 0.1)



# subnetwork cluster
subnet <- induced_subgraph(graph = STRING_net, V(STRING_net)$`display name` %in% hh_results[[1]])
plot(subnet, asp = 0, vertex.size = 2, edge.arrow.size = 0.3, vertex.label.dist = 1)

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
#clust1_net <- simplify(clust1_net)
plot(clust1_net, asp = 0, vertex.label = NA, edge.arrow.size = 0.3, )


# subset edges with experimental evidence
#clust1_net_filt <- delete_edges(clust1_net, E(clust1_net)[!is.na(E(clust1_net)$`stringdb::experiments`) & E(clust1_net)$`stringdb::experiments` >= 0.4])
#plot(clust1_net_filt, asp = 0, vertex.label = NA, vertex.size = 2, edge.arrow.size = 0.3)

# delete the ENSP IDs for later, bites you in the ass with MCODE
subnet <- delete_vertex_attr(subnet, "shared name")
clust1_net <- delete_vertex_attr(clust1_net, "shared name")

write_graph(clust1_net, "WGCNA_subset/results/hhnet_cluster1_netNeighs.graphml", format = "graphml")
write_graph(subnet, "WGCNA_subset/results/hhnet_cluster1_net.graphml", format = "graphml")

clust1_net <- read_graph("WGCNA_subset/results/hhnet_cluster1_netNeighs.graphml", format = "graphml")
subnet <- read_graph("WGCNA_subset/results/hhnet_cluster1_net.graphml", format = "graphml")

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



targets <- read.csv("../Druggability_analysis/data_general/target_all_dbs.csv")
targets <- unique(targets$drugBank_target)

df_subnet <- df_subnet[df_subnet$display.name %in% targets, ]
df_subnetNeighs <- df_subnetNeighs[df_subnetNeighs$display.name %in% targets, ]


