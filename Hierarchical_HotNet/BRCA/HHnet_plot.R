library(tidyverse)
library(igraph)

STRING_net <- read_graph("../../../../OneDrive - RMIT University/PhD/large_git_files/HHnet/STRING_hsa_physical_network.graphml", format = "graphml")
STRING_net <- as.undirected(STRING_net)

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


write_graph(clust1_net, "results/hhnet_cluster1_netNeighs.graphml", format = "graphml")
write_graph(subnet, "results/hhnet_cluster1_net.graphml", format = "graphml")


df_subnet <- data.frame(display.name = V(subnet)$`display name`,
                        degree = degree(subnet),
                        betweenness = betweenness(subnet),
                        source = V(subnet)$color)
df_subnet$source <- ifelse(df_subnet$source == "tomato", "subnet", "STRING")
df_subnet <- df_subnet[order(-df_subnet$degree), ]


df_subnetNeighs <- data.frame(display.name = V(clust1_net)$`display name`,
                              degree = degree(clust1_net),
                              betweenness = betweenness(clust1_net),
                              source = V(clust1_net)$color)
df_subnetNeighs$source <- ifelse(df_subnetNeighs$source == "tomato", "subnet", "STRING")
df_subnetNeighs <- df_subnetNeighs[order(-df_subnetNeighs$degree), ]


# add druggability data
Fpocket_scores <- read.csv("../Druggability_analysis/Fpocket/results_2024.05/fpocket_druggability.csv")

library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

IDs_converted <- getBM(attributes = c( "uniprotswissprot", "external_gene_name", "description"), 
                       filters = "uniprotswissprot", 
                       values = Fpocket_scores$uniprot_id, 
                       mart = ensembl)

unmapped <- IDs_converted[IDs_converted$external_gene_name == "", ]
unrecognised <- Fpocket_scores[!Fpocket_scores$uniprot_id %in% IDs_converted$uniprotswissprot, ]

IDs_converted <- IDs_converted[IDs_converted$external_gene_name != "", ]

novel_transcripts <- unmapped[grep("novel transcript", unmapped$description), ]
novel_proteins <- unmapped[grep("novel protein", unmapped$description), ]
pseudogene <- unmapped[grep("pseudogene", unmapped$description), ]





