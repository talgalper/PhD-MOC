
library(tidyverse)

# load in BRCA data
load("../BRCA_pipe/latest_run/RData/STRING_PPI_FULL.RData")
load("../BRCA_pipe/latest_run/RData/DE_data_geneSymbol.RData")

final_df <- final_df[, -3]

# create index and indexed edge list
gene_index <- as.data.frame(unique(c(final_df$node_1, final_df$node_2)))
gene_index$row_num <- seq.int(nrow(gene_index))
colnames(gene_index) <- c("ensembl_id", "index")
gene_index <- gene_index[c("index", "ensembl_id")]

edge_list_index <- data.frame(from = match(final_df$node_1, gene_index$ensembl_id),
                              to = match(final_df$node_2, gene_index$ensembl_id))

# create score file
DE_data <- subset(DE_data_geneSymbol, select = c("external_gene_name", "logFC"))
DE_data_abs <- subset(DE_data_geneSymbol, select = c("external_gene_name", "logFC_abs"))


write_tsv(gene_index, "data/gene_index.tsv", col_names = F)
write_tsv(edge_list_index, "data/edge_list_index.tsv", col_names = F)
write_tsv(final_df, "data/STRING_network.tsv", col_names = F)
write_tsv(DE_data, "data/logFC_scores.tsv", col_names = F)
write_tsv(DE_data_abs, "data/logFC_scores_abs.tsv", col_names = F)

