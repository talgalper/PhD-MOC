library(tidyverse)
library(igraph)
library(data.table)

STRING_net <- fread("STRING_data/STRING_physical_ENSG.csv")
STRING_net <- STRING_net[!duplicated(t(apply(STRING_net, 1, sort))), ] # get rid of doubled up edges i.e. (a,b) (b,a) same edge weight
STRING_net <- graph_from_data_frame(STRING_net, directed = F)


hh_results <- read_lines("BRCA_CNV/results/clusters_STRING_BRCA_CNV_freq.tsv", skip = 7)
hh_results <- str_split(hh_results, pattern = "\t")

top <- unlist(hh_results[1:10])

#sum(hh_results[[1]] %in% V(STRING_net)$`display name`) 

# colour cluster members
V(STRING_net)$color <- ifelse(V(STRING_net)$name %in% hh_results[[1]], "tomato", "white")
V(STRING_net)$color <- ifelse(V(STRING_net)$name %in% hh_results[[2]], "springgreen", V(STRING_net)$color)
V(STRING_net)$color <- ifelse(V(STRING_net)$name %in% hh_results[[3]], "royalblue", V(STRING_net)$color)
V(STRING_net)$color <- ifelse(V(STRING_net)$name %in% hh_results[[4]], "maroon1", V(STRING_net)$color)
V(STRING_net)$color <- ifelse(V(STRING_net)$name %in% hh_results[[5]], "gold", V(STRING_net)$color)
V(STRING_net)$color <- ifelse(V(STRING_net)$name %in% hh_results[[6]], "orchid", V(STRING_net)$color)
V(STRING_net)$color <- ifelse(V(STRING_net)$name %in% hh_results[[7]], "cyan", V(STRING_net)$color)
V(STRING_net)$color <- ifelse(V(STRING_net)$name %in% hh_results[[8]], "yellowgreen", V(STRING_net)$color)
V(STRING_net)$color <- ifelse(V(STRING_net)$name %in% hh_results[[9]], "mediumseagreen", V(STRING_net)$color)
V(STRING_net)$color <- ifelse(V(STRING_net)$name %in% hh_results[[10]], "saddlebrown", V(STRING_net)$color)



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
subnet <- induced_subgraph(graph = STRING_net, V(STRING_net)$name %in% top)
set.seed(1234)
plot.igraph(subnet, asp = 0, vertex.label = NA, vertex.size = 2, edge.arrow.size = 0.3)

# add direct neighbors of cluster nodes
neighbs <- neighborhood(
  STRING_net,
  order = 1,
  nodes = top
)
neighbs <- flatten(neighbs)
neighbs <- names(neighbs)
neighbs <- unique(neighbs)

clust1_net <- induced_subgraph(STRING_net, V(STRING_net)[V(STRING_net)$name %in% neighbs])
V(clust1_net)$size <- scales::rescale(degree(clust1_net), c(1, 7))
#clust1_net <- simplify(clust1_net) # pretty sure gets rid of all extra data from STRING/Cytoscape
set.seed(1234)
plot.igraph(clust1_net, asp = 0, vertex.label = NA, edge.arrow.size = 0.3, edge.curved = F)


# subset edges with experimental evidence
#clust1_net_filt <- delete_edges(clust1_net, E(clust1_net)[!is.na(E(clust1_net)$`stringdb::experiments`) & E(clust1_net)$`stringdb::experiments` >= 0.4])
#plot(clust1_net_filt, asp = 0, vertex.label = NA, vertex.size = 2, edge.arrow.size = 0.3)


write_graph(clust1_net, "BRCA_CNV/results/hhnet_cluster1_netNeighs.graphml", format = "graphml")
write_graph(subnet, "BRCA_CNV/results/hhnet_cluster1_net.graphml", format = "graphml")

clust1_net <- read_graph("BRCA_CNV/results/hhnet_cluster1_netNeighs.graphml", format = "graphml")
subnet <- read_graph("BRCA_CNV/results/hhnet_cluster1_net.graphml", format = "graphml")


df_subnet <- data.frame(ENSG = V(subnet)$name,
                        degree = degree(subnet),
                        betweenness = betweenness(subnet),
                        closeness = closeness(subnet),
                        eigen_centrality = eigen_centrality(subnet)$vector,
                        cluster = V(subnet)$color)
df_subnet$source <- ifelse(df_subnet$cluster != "white", "subnet", "STRING")


library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

ensembl_converted <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description", "gene_biotype"), 
                           filters = "ensembl_gene_id", 
                           values = df_subnet$ENSG, 
                           mart = ensembl)
ensembl_converted$description <- gsub("\\[.*?\\]", "", ensembl_converted$description)

unmapped <- ensembl_converted[ensembl_converted$external_gene_name == "", ]
unrecognised <- df_subnet[!df_subnet$ENSG %in% ensembl_converted$ensembl_gene_id, ]

df_subnet <- merge.data.table(ensembl_converted, df_subnet, by.x = "ensembl_gene_id", by.y = "ENSG", all.y = T)
df_subnet <- df_subnet[order(-df_subnet$degree), ]
rownames(df_subnet) <- NULL



df_subnetNeighs <- data.frame(ENSG = V(clust1_net)$name,
                              degree = degree(clust1_net),
                              betweenness = betweenness(clust1_net),
                              closeness = closeness(clust1_net),
                              eigen_centrality = eigen_centrality(clust1_net)$vector,
                              cluster = V(clust1_net)$color)
df_subnetNeighs$source <- ifelse(df_subnetNeighs$cluster != "white", "subnet", "STRING")

ensembl_converted <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description", "gene_biotype"), 
                           filters = "ensembl_gene_id", 
                           values = df_subnetNeighs$ENSG, 
                           mart = ensembl)
ensembl_converted$description <- gsub("\\[.*?\\]", "", ensembl_converted$description)

unmapped <- ensembl_converted[ensembl_converted$external_gene_name == "", ]
unrecognised <- df_subnetNeighs[!df_subnetNeighs$ENSG %in% ensembl_converted$ensembl_gene_id, ]
df_subnetNeighs <- merge.data.table(ensembl_converted, df_subnetNeighs, by.x = "ensembl_gene_id", by.y = "ENSG", all.y = T)

df_subnetNeighs <- df_subnetNeighs[order(-df_subnetNeighs$degree), ]
rownames(df_subnetNeighs) <- NULL

fwrite(df_subnet, "BRCA_CNV/results/df_subnet.csv")
fwrite(df_subnetNeighs, "BRCA_CNV/results/df_subnetNeighs.csv")



# network metrics of known BRCA targets
df_subnet <- fread("BRCA_CNV/results/df_subnet.csv")
df_subnetNeighs <- fread("BRCA_CNV/results/df_subnetNeighs.csv")

targets <- read.csv("../Druggability_analysis/data_general/target_all_dbs.csv")
targets <- unique(targets$ensembl_gene_id)

df_subnet_targets <- df_subnet[df_subnet$ensembl_gene_id %in% targets, ]
df_subnetNeighs_targets <- df_subnetNeighs[df_subnetNeighs$ensembl_gene_id %in% targets, ]




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



library(clusterProfiler)
library(progress)
pb <- progress_bar$new(
  format = "  Performing GO Analysis [:bar] :percent eta: :eta",
  total = 10
)

cluster_GO <- list()
clust_no <- 1
for (cluster in hh_results[1:18]) {
  GO <- enrichGO(cluster, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")
  cluster_GO[[clust_no]] <- GO
  clust_no <- clust_no + 1
  rm(GO, cluster)
  pb$tick()
}
rm(clust_no, pb)

GO_formatted <- data.frame()
clust_no <- 1
for (cluster in cluster_GO) {
  result <- cluster@result
  result_top <- head(result, 10)
  result_top$cluster <- rep(clust_no, nrow(result_top))
  
  GO_formatted <- rbind(GO_formatted, result_top)
  rm(result_top, result)
  clust_no <- clust_no + 1
}
rm(clust_no)

GO_formatted$color <- ifelse(GO_formatted$cluster == 1, "tomato", "")
GO_formatted$color <- ifelse(GO_formatted$cluster == 2, "springgreen", GO_formatted$color)
GO_formatted$color <- ifelse(GO_formatted$cluster == 3, "royalblue", GO_formatted$color)
GO_formatted$color <- ifelse(GO_formatted$cluster == 4, "maroon1", GO_formatted$color)
GO_formatted$color <- ifelse(GO_formatted$cluster == 5, "gold", GO_formatted$color)
GO_formatted$color <- ifelse(GO_formatted$cluster == 6, "orchid", GO_formatted$color)
GO_formatted$color <- ifelse(GO_formatted$cluster == 7, "cyan", GO_formatted$color)
GO_formatted$color <- ifelse(GO_formatted$cluster == 8, "yellowgreen", GO_formatted$color)
GO_formatted$color <- ifelse(GO_formatted$cluster == 9, "mediumseagreen", GO_formatted$color)
GO_formatted$color <- ifelse(GO_formatted$cluster == 10, "saddlebrown", GO_formatted$color)

save(cluster_GO, file = "~/OneDrive - RMIT University/PhD/large_git_files/HHnet/HHnet_cluster_GO.RData")
load("~/OneDrive - RMIT University/PhD/large_git_files/HHnet/HHnet_cluster_GO.RData")
fwrite(GO_formatted, "BRCA/STN_filt/results/subnet_cluster_GO.csv")



subnet_GO <- enrichGO(df_subnet$ensembl_gene_id, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")
subnet_GO <- simplify(subnet_GO)
subnet_GO <- subnet_GO@result

subnetNeighs_GO <- enrichGO(df_subnetNeighs$ensembl_gene_id, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")
subnetNeighs_GO <- simplify(subnetNeighs_GO)
subnetNeighs_GO <- subnetNeighs_GO@result

fwrite(subnet_GO, "BRCA/STN_filt/results/subnet_GO.csv")
fwrite(subnetNeighs_GO, "BRCA/STN_filt/results/subnetNeighs_GO.csv")





BRCA_markers <- c("ERBB2", "TP53", "BRCA1", "BRCA2", "ESR1", "PGR", "PIK3CA", "PTEN",
                  "MYC", "GATA3", "CDH1", "AKT1", "CCND1", "FGFR1", "MUC1", "EGFR",
                  "NOTCH1", "KRAS", "BCL2", "MAPK1")

temp <- df_subnetNeighs[df_subnetNeighs$external_gene_name %in% BRCA_markers, ]



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












# add direct neighbors of cluster nodes
neighbs <- neighborhood(
  STRING_net,
  order = 1,
  nodes = hh_results[[3]]
)
neighbs <- flatten(neighbs)
neighbs <- names(neighbs)
neighbs <- unique(neighbs)

clust1_net <- induced_subgraph(STRING_net, V(STRING_net)[V(STRING_net)$name %in% neighbs])
V(clust1_net)$size <- scales::rescale(degree(clust1_net), c(1, 7))
#clust1_net <- simplify(clust1_net) # pretty sure gets rid of all extra data from STRING/Cytoscape
set.seed(1234)
plot.igraph(clust1_net, asp = 0, vertex.label = NA, edge.arrow.size = 0.3, edge.curved = F)



df_subnetNeighs <- data.frame(ENSG = V(clust1_net)$name,
                              degree = degree(clust1_net),
                              betweenness = betweenness(clust1_net),
                              closeness = closeness(clust1_net),
                              eigen_centrality = eigen_centrality(clust1_net)$vector,
                              cluster = V(clust1_net)$color)
df_subnetNeighs$source <- ifelse(df_subnetNeighs$cluster != "white", "subnet", "STRING")

ensembl_converted <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description", "gene_biotype"), 
                           filters = "ensembl_gene_id", 
                           values = df_subnetNeighs$ENSG, 
                           mart = ensembl)
ensembl_converted$description <- gsub("\\[.*?\\]", "", ensembl_converted$description)

unmapped <- ensembl_converted[ensembl_converted$external_gene_name == "", ]
unrecognised <- df_subnetNeighs[!df_subnetNeighs$ENSG %in% ensembl_converted$ensembl_gene_id, ]
df_subnetNeighs <- merge.data.table(ensembl_converted, df_subnetNeighs, by.x = "ensembl_gene_id", by.y = "ENSG", all.y = T)

df_subnetNeighs <- df_subnetNeighs[order(-df_subnetNeighs$degree), ]
rownames(df_subnetNeighs) <- NULL




targets <- read.csv("../Druggability_analysis/data_general/target_all_dbs.csv")
targets <- unique(targets$ensembl_gene_id)

df_subnet_targets <- df_subnet[df_subnet$ensembl_gene_id %in% targets, ]
df_subnetNeighs_targets <- df_subnetNeighs[df_subnetNeighs$ensembl_gene_id %in% targets, ]




