library(PCSF)
library(tidyverse)
library(edgeR)

# Load data
load("RData/LumA/DE_data.RData")
load("RData/LumB/DE_data.RData")
load("RData/Her2/DE_data.RData")
load("RData/basal/DE_data.RData")
load("RData/TCGA_normal.RData")

GTEx_data <- read.table("../BRCA_pipe/gene_reads_2017-06-05_v8_breast_mammary_tissue.gct", skip = 2)
colnames(GTEx_data) <- GTEx_data[1, ]
GTEx_data <- GTEx_data[-1, -1]
rownames(GTEx_data) <- NULL

# opt for having gene Ensembl IDs instead of gene names as rownames (same as TCGA)
GTEx_ENS <- column_to_rownames(GTEx_data, "Name")
rownames(GTEx_ENS) <- gsub("\\.\\d+", "", rownames(GTEx_ENS))
GTEx_ENS <- GTEx_ENS[ , -1]
rownames <- rownames(GTEx_ENS)
GTEx_ENS <- as.data.frame(sapply(GTEx_ENS, as.numeric))
rownames(GTEx_ENS) <- rownames
rm(rownames, GTEx_data)
GTEx_ENS[] <- lapply(GTEx_ENS, function(x){as.integer(x)})

# combine all tumour samples
all_subtypes <- cbind(LumA_unstranded, LumB_unstranded, Her2_unstranded, Basal_unstranded)

# clean env
rm(normal_unstranded, LumA_unstranded, LumB_unstranded, Her2_unstranded, Basal_unstranded)
collectGarbage()

# read in functions from "../BRCA_pipe/Rscript/DE_functions.R"
counts_filt <- filter_low_expr(disease_data = all_subtypes, 
                                control_data = GTEx_ENS)
hist(log(as.matrix(counts_filt$counts_filt)))

DE_results <- DE_analysis(counts_matrix = DE_counts_filt$counts_filt,
                          sample_info = DE_counts_filt$sample_info)

# load in DE results from WGCNA where we already did this analysis for BRCA
load("../WGCNA/BRCA/RData/DE_subset/dif_exp.RData")

DE_data <- subset(dif_exp, select = c("gene_id", "logFC"))
DE_data$logFC_abs <- abs(DE_data$logFC) # get absolute values

write.table(DE_data$gene_id, "latest_run/intermediate/dif_exp_genes.txt", col.names = F, row.names = F, quote = F)


# get interaction data from STRING
string_edge_data <- read.table("latest_run/intermediate/STRING network (physical) default edge.csv", header = T, sep = ",", stringsAsFactors = F)
ppi_list <- subset(string_edge_data, select = c("name", "stringdb..score"))
ppi_list <- ppi_list %>% 
  separate(name, sep = " ", into = c("node_1", "del", "node_2"))
ppi_list <- subset(ppi_list, select = c("node_1", "node_2", "stringdb..score"))
ppi_list$node_1 <- gsub(".*.\\.", "", ppi_list$node_1)
ppi_list$node_2 <- gsub(".*.\\.", "", ppi_list$node_2)

string_node_data <- read.table("latest_run/intermediate/STRING network (physical) default node.csv", header = T, sep = ",", stringsAsFactors = F)
node_list <- subset(string_node_data, select = c("name", "query.term"))
node_list$name <- gsub(".*.\\.", "", node_list$name)
ppi_list$original_order <- seq_len(nrow(ppi_list))
merged_df <- merge(ppi_list, node_list, by.x = "node_1", by.y = "name", all.x = TRUE)
merged_df <- merge(merged_df, node_list, by.x = "node_2", by.y = "name", all.x = TRUE)
merged_df <- merged_df[order(merged_df$original_order), ]

final_df <- merged_df[, c("query.term.x", "query.term.y", "stringdb..score")]
colnames(final_df) <- c("node_1", "node_2", "score")

save(final_df, DE_data, file = "latest_run/RData/PCSF_input.RData")



# master human STRING network
#STRING <- read.table("latest_run/intermediate/9606.protein.physical.links.v12.0.txt", header = T)
#STRING_full <- read.table("latest_run/intermediate/9606.protein.physical.links.full.v12.0.txt", header = T)
#STRING_info <- read.table("latest_run/intermediate/9606.protein.info.v12.0.txt", header = T, sep = "\t")
#STRING_aliases <- read.table("latest_run/intermediate/9606.protein.aliases.v12.0.txt", header = T, "\t", fill = T)
#STRING_aliases <- STRING_aliases[STRING_aliases$source == "Ensembl_gene", ]
#
## Replace with ensembl gene
#STRING_updated <- STRING_full
#STRING_updated <- STRING_updated[STRING_updated$experiments != 0, ]
#STRING_updated <- STRING_updated[STRING_updated$combined_score >= 200, ]
#
#
#STRING_updated$protein1_gene <- ifelse(is.na(STRING_aliases$alias[match(STRING_updated$protein1, STRING_aliases$string_protein_id)]),
#                                       STRING_updated$protein1,
#                                       STRING_aliases$alias[match(STRING_updated$protein1, STRING_aliases$string_protein_id)])
#
#STRING_updated$protein2_gene <- ifelse(is.na(STRING_aliases$alias[match(STRING_updated$protein2, STRING_aliases$string_protein_id)]),
#                                       STRING_updated$protein2,
#                                       STRING_aliases$alias[match(STRING_updated$protein2, STRING_aliases$string_protein_id)])
#
#STRING_updated <- STRING_updated[, c("protein1_gene", "protein2_gene", "combined_score")]
#
#STRING_updated$protein1_gene <- gsub(".*.\\.", "", STRING_updated$protein1_gene)
#STRING_updated$protein2_gene <- gsub(".*.\\.", "", STRING_updated$protein2_gene)


# try PCSF with logFC scores from WGCNA results
load("../WGCNA/BRCA/RData/all_default/venn_data.RData")
common_genes <- Reduce(intersect, list(DE_genes, tumour_associated, top_kwithin, top_gene_membership))
DE_data_subset <- DE_data[DE_data$gene_id %in% common_genes, ]



load("latest_run/RData/PCSF_input.RData")

# set seed for reproducibility 
set.seed(1234)
# construct interactome
ppi <- construct_interactome(final_df)
# set terminals
terminals <- setNames(as.numeric(DE_data_subset$logFC_abs), DE_data_subset$gene_id)

# run PCSF with random noise
start_time <- Sys.time()
subnet <- PCSF_rand(ppi, terminals, n = 50, r = 0.1, b = 1, w = 2, mu = 0.0005)
elapsed_time <- Sys.time() - start_time
print(elapsed_time)

save(subnet, file = "latest_run/RData/PCSF_subnet.RData")

# convert to gene symbols
library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene_names <- V(subnet)$name

ensembl_converted <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), 
                           filters = "ensembl_gene_id", 
                           values = gene_names, 
                           mart = ensembl)

V(subnet)$name <- ensembl_converted$external_gene_name

plot.PCSF(subnet, node_label_cex = 15)

# enrich network using enrichR
PCSF_enriched <- enrichment_analysis(subnet)
plot.PCSFe(PCSF_enriched$subnet, edge_width = 8, node_size = 30, node_label_cex = 1)


load("RData/LumA/filterByExp/GTEx/PCSF_subnet.RData")
load("RData/LumB/filterByExp/GTEx/PCSF_subnet.RData")
load("RData/Her2/filterByExp/GTEx/PCSF_subnet.RData")
load("RData/basal/filterByExp/GTEx/PCSF_subnet.RData")


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

rownames(df) <- 1:nrow(df)

df <- df[order(-df$degree_centrality), ]
rownames(df) <- NULL


ensembl_converted <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), 
                           filters = "ensembl_gene_id", 
                           values = df$gene_id, 
                           mart = ensembl)

df <- merge(df, ensembl_converted, by.x = "gene_id", by.y = "ensembl_gene_id", all.x = T)
