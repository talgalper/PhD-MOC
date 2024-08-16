library(PCSF)
library(tidyverse)
library(edgeR)


# load in DE results from WGCNA where we already did this analysis for BRCA
load("../WGCNA/MOC/RData/MOC_dif_exp.RData")
DE_data <- subset(dif_exp, select = c("gene_id", "logFC"))
DE_data$logFC_abs <- abs(DE_data$logFC) # get absolute values

# convert to gene symbols
library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

ensembl_converted <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description"), 
                           filters = "ensembl_gene_id", 
                           values = DE_data$gene_id, 
                           mart = ensembl)

unmapped <- ensembl_converted[ensembl_converted$external_gene_name == "", ]
unrecognised <- DE_data[!DE_data$gene_id %in% ensembl_converted$ensembl_gene_id, ]

ensembl_converted <- ensembl_converted[ensembl_converted$external_gene_name != "", ]

novel_transcripts <- unmapped[grep("novel transcript", unmapped$description), ]
novel_proteins <- unmapped[grep("novel protein", unmapped$description), ]
pseudogene <- unmapped[grep("pseudogene", unmapped$description), ]

DE_data_geneSymbol <- merge(DE_data, ensembl_converted, by.x = "gene_id", by.y = "ensembl_gene_id")
DE_data_geneSymbol <- subset(DE_data_geneSymbol, select = c("external_gene_name", "logFC", "logFC_abs"))

save(DE_data_geneSymbol, unmapped, unrecognised, file = "latest_run/MOC/DE_data_geneSymbol.RData")

# get interaction data from STRING
string_edge_data <- read.table("latest_run/intermediate/STRING network (physical) default edge - FULL.csv", header = T, sep = ",", stringsAsFactors = F)
ppi_list <- subset(string_edge_data, select = c("name", "stringdb..score"))
ppi_list <- ppi_list %>% 
  separate(name, sep = " ", into = c("node_1", "del", "node_2"))
ppi_list <- subset(ppi_list, select = c("node_1", "node_2", "stringdb..score"))
ppi_list$node_1 <- gsub(".*.\\.", "", ppi_list$node_1)
ppi_list$node_2 <- gsub(".*.\\.", "", ppi_list$node_2)

string_node_data <- read.table("latest_run/intermediate/STRING network (physical) default node - FULL.csv", header = T, sep = ",", stringsAsFactors = F)
node_list <- subset(string_node_data, select = c("name", "display.name"))
node_list$name <- gsub(".*.\\.", "", node_list$name)
ppi_list$original_order <- seq_len(nrow(ppi_list))
merged_df <- merge(ppi_list, node_list, by.x = "node_1", by.y = "name", all.x = TRUE)
merged_df <- merge(merged_df, node_list, by.x = "node_2", by.y = "name", all.x = TRUE)
merged_df <- merged_df[order(merged_df$original_order), ]

final_df <- merged_df[, c("display.name.x", "display.name.y", "stringdb..score")]
colnames(final_df) <- c("node_1", "node_2", "score")

save(final_df, file = "latest_run/RData/STRING_PPI_FULL.RData")
load("latest_run/RData/STRING_PPI_FULL.RData")


# master human STRING network
#STRING <- read.table("latest_run/intermediate/9606.protein.physical.links.v12.0.txt", header = T)
#STRING_full <- read.table("latest_run/intermediate/9606.protein.physical.links.full.v12.0.txt", header = T)
#STRING_info <- read.table("latest_run/intermediate/9606.protein.info.v12.0.txt", header = T, sep = "\t")
#STRING_aliases <- read.table("../../../../OneDrive - RMIT University/PhD/large_git_files/9606.protein.aliases.v12.0.txt", header = T, "\t", fill = T)
#STRING_aliases <- STRING_aliases[STRING_aliases$source == "Ensembl_gene", ]
#
## Replace with preffered name
#STRING_updated <- STRING_full
#STRING_updated <- STRING_updated[STRING_updated$experiments != 0, ]
#STRING_updated <- STRING_updated[STRING_updated$combined_score >= 200, ]
#
#STRING_updated$protein1_gene <- STRING_info$preferred_name[match(STRING_updated$protein1, STRING_info$string_protein_id)]
#STRING_updated$protein2_gene <- STRING_info$preferred_name[match(STRING_updated$protein2, STRING_info$string_protein_id)]
#
#STRING_updated <- STRING_updated[, c("protein1_gene", "protein2_gene", "combined_score")]
#STRING_updated <- na.omit(STRING_updated)



# try PCSF with logFC scores from WGCNA results
load("../WGCNA/BRCA/RData/all_default/venn_data.RData")
common_genes <- Reduce(intersect, list(DE_genes, tumour_associated, top_kwithin, top_gene_membership))
DE_data_subset <- DE_data[DE_data$gene_id %in% common_genes, ]



# set seed for reproducibility 
set.seed(1234)
# construct interactome
ppi <- construct_interactome(final_df)
# set terminals
terminals <- setNames(as.numeric(DE_data_geneSymbol$logFC_abs), DE_data_geneSymbol$external_gene_name)

# run PCSF with random noise
start_time <- Sys.time()
subnet <- PCSF_rand(ppi, terminals, n = 50, r = 0.1, b = 1, w = 2, mu = 0.0005)
elapsed_time <- Sys.time() - start_time
print(elapsed_time)

plot.PCSF(subnet, node_label_cex = 15)

save(subnet, file = "latest_run/MOC/PCSF_MOCsubnet_FULL.RData")
load("latest_run/MOC/PCSF_MOCsubnet_FULL.RData")


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


# network enrichment
PCSF_enrich <- enrichment_analysis(subnet)

# Create a data frame with the enrichment results
enrichment_results <- PCSF_enrich$enrichment
enrichment_table <- data.frame(
  Cluster = enrichment_results$Cluster,
  Term = enrichment_results$Term,
  PValue = enrichment_results$P.value,
  Adjusted_Pvalue = enrichment_results$Adjusted.P.value,
  Genes = enrichment_results$Genes)
