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
df$closeness <- closeness(subnet)
df$prize <- V(subnet)$prize
df$type <- V(subnet)$type

rownames(df) <- 1:nrow(df)

df <- df[order(-df$degree_centrality), ]
rownames(df) <- NULL



# add gene symbols
library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

ensembl_converted <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description", "gene_biotype"), 
                           filters = "ensembl_gene_id", 
                           values = df$gene_id, 
                           mart = ensembl)
ensembl_converted$description <- gsub("\\[.*?\\]", "", ensembl_converted$description)

unmapped <- ensembl_converted[ensembl_converted$external_gene_name == "", ]
unrecognised <- df[!df$gene_id %in% ensembl_converted$ensembl_gene_id, ]

df <- merge.data.table(ensembl_converted, df, by.x = "ensembl_gene_id", by.y = "gene_id", all.y = T)
df <- df[order(-df$degree_centrality), ]
rownames(df) <- NULL

save(df, file = "latest_run/RData/STN_filt/PCSF_results.RData")
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



