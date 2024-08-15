
library(tidyverse)

# load in data
load("../../Documents/GitHub/PhD-MOC/BRCA_pipe/latest_run/RData/STRING_PPI_FULL.RData")
load("../../Documents/GitHub/PhD-MOC/WGCNA/MOC/RData/MOC_dif_exp.RData")

library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

ensembl_converted <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description"), 
                           filters = "ensembl_gene_id", 
                           values = dif_exp$gene_id, 
                           mart = ensembl)

unmapped <- ensembl_converted[ensembl_converted$external_gene_name == "", ]
unrecognised <- dif_exp[!dif_exp$gene_id %in% ensembl_converted$ensembl_gene_id, ]

ensembl_converted <- ensembl_converted[ensembl_converted$external_gene_name != "", ]

DE_data_geneSymbol <- merge(dif_exp, ensembl_converted, by.x = "gene_id", by.y = "ensembl_gene_id")
DE_data_geneSymbol <- subset(DE_data_geneSymbol, select = c("external_gene_name", "logFC"))

DE_data_geneSymbol$logFC_abs <- abs(DE_data_geneSymbol$logFC) # get absolute values


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


write_tsv(gene_index, "MOC/data/gene_index.tsv", col_names = F)
write_tsv(edge_list_index, "MOC/data/edge_list_index.tsv", col_names = F)
write_tsv(final_df, "MOC/data/STRING_network.tsv", col_names = F)
write_tsv(DE_data, "MOC/data/logFC_scores.tsv", col_names = F)
write_tsv(DE_data_abs, "MOC/data/logFC_scores_abs.tsv", col_names = F)

