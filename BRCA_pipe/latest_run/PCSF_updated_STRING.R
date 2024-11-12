library(data.table)
library(PCSF)


# load in DE results from WGCNA where we already did this analysis for BRCA
load("../WGCNA/BRCA/RData/DE_subset/dif_exp.RData")
DE_data <- subset(dif_exp, select = c("gene_id", "logFC"))
DE_data$logFC_abs <- abs(DE_data$logFC) # get absolute values

STRING_edge <- fread("latest_run/intermediate/STRING_physical_ENSG.csv")

# set seed for reproducibility 
set.seed(1234)
# construct interactome
ppi <- construct_interactome(STRING_edge)
# set terminals
terminals <- setNames(as.numeric(DE_data$logFC_abs), DE_data$gene_id)

# run PCSF with random noise
start_time <- Sys.time()
subnet <- PCSF_rand(ppi, terminals, n = 50, r = 0.1, b = 1, w = 2, mu = 0.0005)
elapsed_time <- Sys.time() - start_time
print(elapsed_time)

plot.PCSF(subnet, node_label_cex = 15)

save(subnet, file = "latest_run/RData/PCSF_updated_STRING.RData")
load("latest_run/RData/PCSF_updated_STRING.RData")



# extract cluster data
clust <- components(subnet)
df <- data.frame(gene_id = names(clust$membership), cluster = factor(clust$membership))
betweenness <- betweenness(subnet) 
centrality <- degree(subnet) 
df$betweenness <- betweenness[as.character(df$gene_id)]
df$degree_centrality <- centrality[as.character(df$gene_id)]
df$betweenness <- as.integer(df$betweenness)
df$degree_centrality <- as.integer(df$degree_centrality)
df$prize <- V(subnet)$prize
df$type <- V(subnet)$type

rownames(df) <- 1:nrow(df)

df <- df[order(-df$degree_centrality), ]
rownames(df) <- NULL
