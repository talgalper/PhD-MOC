library(data.table)
library(PCSF)


# load in DE results from WGCNA where we already did this analysis for BRCA
load("latest_run/RData/STN_filt/dif_exp.RData")
DE_data <- subset(dif_exp, select = c("gene_id", "logFC"))
DE_data$logFC_abs <- abs(DE_data$logFC) # get absolute values

STRING_edge <- fread("latest_run/intermediate/STRING_physical_ENSG.csv", data.table = F)
STRING_edge <- STRING_edge[!duplicated(t(apply(STRING_edge, 1, sort))), ]


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

save(subnet, file = "latest_run/RData/STN_filt/PCSF_updated_STRING.RData")
load("latest_run/RData/STN_filt/PCSF_updated_STRING.RData")



# extract cluster data
clust <- components(subnet)
df <- data.frame(gene_id = names(clust$membership), cluster = factor(clust$membership))
betweenness <- betweenness(subnet) 
centrality <- degree(subnet) 
df$betweenness <- betweenness[as.character(df$gene_id)]
df$degree_centrality <- centrality[as.character(df$gene_id)]
df$betweenness <- as.integer(df$betweenness)
df$degree_centrality <- as.integer(df$degree_centrality)
df$eigen_centrality <- eigen_centrality(subnet)$vector
df$page_rank <- page_rank(subnet)$vector
df$closeness <- closeness(subnet)
df$prize <- V(subnet)$prize
df$type <- V(subnet)$type

df <- df[order(-df$degree_centrality), ]
rownames(df) <- NULL



# add gene symbols - read in function from "../MOC_pipe/R/functions.R"
df <- id_annot(data = df,
               input_type = "ensembl_gene_id", 
               convert_to = c("external_gene_name", "description", "gene_biotype"))

df <- df[order(-df$degree_centrality), ]
rownames(df) <- NULL

save(df, file = "latest_run/RData/STN_filt/PCSF_results_pageRank.RData")
load("latest_run/RData/STN_filt/PCSF_results.RData")

library(clusterProfiler)
PCSF_GO <- enrichGO(df$ensembl_gene_id, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")
PCSF_GO <- simplify(PCSF_GO)
PCSF_GO <- PCSF_GO@result

save(PCSF_GO, file = "latest_run/RData/STN_filt/PCSF_GO.RData")

targets <- read.csv("../Druggability_analysis/data_general/target_all_dbs.csv")
targets <- targets[, c(2,4)]
targets <- targets[!duplicated(targets$ensembl_gene_id), ]

temp <- df[df$ensembl_gene_id %in% targets$ensembl_gene_id, ]



