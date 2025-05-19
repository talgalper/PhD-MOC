library(igraph)
library(data.table)

STRING_net <- fread("STRING_data/STRING_physical_ENSG.csv")
STRING_net <- STRING_net[!duplicated(t(apply(STRING_net, 1, sort))), ] # get rid of doubled up edges i.e. (a,b) (b,a) same edge weight
STRING_net <- graph_from_data_frame(STRING_net, directed = F)


hh_results <- readr::read_lines("MOC/results/clusters_STRING_MOC_logFC.tsv", skip = 7)
hh_results <- stringr::str_split(hh_results, pattern = "\t")
table(lengths(hh_results))

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

subnet_neighs <- induced_subgraph(STRING_net, V(STRING_net)[V(STRING_net)$name %in% neighbs])
V(subnet_neighs)$size <- scales::rescale(degree(subnet_neighs), c(1, 7))
#clust1_net <- simplify(clust1_net) # pretty sure gets rid of all extra data from STRING/Cytoscape
set.seed(1234)
plot.igraph(subnet_neighs, asp = 0, vertex.label = NA, edge.arrow.size = 0.3, edge.curved = F)


# subset edges with experimental evidence
#clust1_net_filt <- delete_edges(clust1_net, E(clust1_net)[!is.na(E(clust1_net)$`stringdb::experiments`) & E(clust1_net)$`stringdb::experiments` >= 0.4])
#plot(clust1_net_filt, asp = 0, vertex.label = NA, vertex.size = 2, edge.arrow.size = 0.3)


write_graph(subnet_neighs, "MOC/results/hhnet_netNeighs.graphml", format = "graphml")
write_graph(subnet, "MOC/results/hhnet_subnet.graphml", format = "graphml")

subnet_neighs <- read_graph("MOC/results/hhnet_netNeighs.graphml", format = "graphml")
subnet <- read_graph("MOC/results/hhnet_subnet.graphml", format = "graphml")


df_subnet <- data.frame(ENSG = V(subnet)$name,
                        degree = degree(subnet),
                        betweenness = betweenness(subnet),
                        closeness = closeness(subnet),
                        eigen_centrality = eigen_centrality(subnet)$vector,
                        page_rank = page_rank(subnet)$vector,
                        cluster = V(subnet)$color)
df_subnet$source <- ifelse(df_subnet$cluster != "white", "subnet", "STRING")


# convert to gene symbols
library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

ensembl_id_annot <- function(ensembl, data, col_id = 1, input_type, convert_to) {
  library(tidyverse)
  if (col_id == 0) {
    data <- rownames_to_column(data)
    col_id <- 1
  }
  
  ensembl_annot <- getBM(attributes = c(input_type, convert_to), 
                         filters = input_type, 
                         values = data[, col_id], 
                         mart = ensembl)
  ensembl_annot$description <- gsub("\\[.*?\\]", "", ensembl_annot$description)
  
  data_annot <- merge(ensembl_annot, data, by.x = input_type, by.y = colnames(data)[col_id], all.y = T)
}

df_subnet <- ensembl_id_annot(ensembl, df_subnet,
                              input_type = "ensembl_gene_id",
                              convert_to = c("external_gene_name", "description", "gene_biotype"))
df_subnet <- df_subnet[order(-df_subnet$degree), ]
rownames(df_subnet) <- NULL


df_subnetNeighs <- data.frame(ENSG = V(subnet_neighs)$name,
                              degree = degree(subnet_neighs),
                              betweenness = betweenness(subnet_neighs),
                              closeness = closeness(subnet_neighs),
                              eigen_centrality = eigen_centrality(subnet_neighs)$vector,
                              page_rank = page_rank(subnet_neighs)$vector,
                              cluster = V(subnet_neighs)$color)
df_subnetNeighs$source <- ifelse(df_subnetNeighs$cluster != "white", "subnet", "STRING")


df_subnetNeighs <- ensembl_id_annot(ensembl, df_subnetNeighs,
                                    input_type = "ensembl_gene_id",
                                    convert_to = c("external_gene_name", "description", "gene_biotype"))
df_subnetNeighs <- df_subnetNeighs[order(-df_subnetNeighs$degree), ]
rownames(df_subnetNeighs) <- NULL


fwrite(df_subnet, "MOC/results/df_subnet.csv")
fwrite(df_subnetNeighs, "MOC/results/df_subnetNeighs.csv")





library(clusterProfiler)
library(progress)

GO <- enrichGO(MOCvsBEN$external_gene_name, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP", universe = rownames(filtered_MOC_data))
result <- GO@result

load("../MOC_pipe/DE/MOC_vs_BEN/DE_results.RData")
filtered_MOC_data <- DE_results$input_data

pb <- progress_bar$new(
  format = "  Performing GO Analysis [:bar] :percent eta: :eta",
  total = 18
)

cluster_GO <- list()
clust_no <- 1
for (cluster in hh_results[1:10]) {
  GO <- enrichGO(cluster, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP", universe = rownames(filtered_MOC_data))
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

save(cluster_GO, file = "~/OneDrive - RMIT University/PhD/large_git_files/MOC/MOCvsBEN_HHnet_cluster_GO.RData")
fwrite(GO_formatted, "MOC/results/subnetNeighs_cluster_GO.csv")

load("~/OneDrive - RMIT University/PhD/large_git_files/HHnet/HHnet_cluster_GO.RData")



