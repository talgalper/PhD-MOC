library(tidyverse)

load("~/Documents/GitHub/PhD-MOC/BRCA_CNV/RData/BRCA_CNV_summary.RData")
load("../BRCA_pipe/latest_run/RData/STRING_PPI_FULL.RData")

final_df <- final_df[, -3]

# create index and indexed edge list
gene_index <- as.data.frame(unique(c(final_df$node_1, final_df$node_2)))
gene_index$row_num <- seq.int(nrow(gene_index))
colnames(gene_index) <- c("ensembl_id", "index")
gene_index <- gene_index[c("index", "ensembl_id")]

edge_list_index <- data.frame(from = match(final_df$node_1, gene_index$ensembl_id),
                              to = match(final_df$node_2, gene_index$ensembl_id))

# create score file
CNV_data <- subset(cnv_summary_geneSymbol, select = c("external_gene_name", "alterations_percentage"))
CNV_data$alterations_percentage <- CNV_data$alterations_percentage / 100

write_tsv(gene_index, "BRCA_CNV/data/gene_index.tsv", col_names = F)
write_tsv(edge_list_index, "BRCA_CNV/data/edge_list_index.tsv", col_names = F)
write_tsv(final_df, "BRCA_CNV/data/STRING_network.tsv", col_names = F)
write_tsv(CNV_data, "BRCA_CNV/data/CNV_scores.tsv", col_names = F)

