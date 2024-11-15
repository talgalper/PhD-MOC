library(data.table)
library(tidyverse)

# load in BRCA data
STRING_edge <- fread("STRING_data/STRING_physical_ENSG.csv") 
STRING_edge <- STRING_edge[, -3]
STRING_edge <- STRING_edge[!duplicated(t(apply(STRING_edge, 1, sort))), ]


# create index and indexed edge list
gene_index <- as.data.frame(unique(c(STRING_edge$protein1_ENSG, STRING_edge$protein2_ENSG)))
gene_index$row_num <- seq.int(nrow(gene_index))
colnames(gene_index) <- c("ensembl_id", "index")
gene_index <- gene_index[c("index", "ensembl_id")]

edge_list_index <- data.frame(from = match(STRING_edge$protein1_ENSG, gene_index$ensembl_id),
                              to = match(STRING_edge$protein2_ENSG, gene_index$ensembl_id))

# create score file
load("../BRCA_pipe/latest_run/RData/STN_filt/dif_exp.RData")
DE_data <- subset(dif_exp, select = c("gene_id", "logFC"))
DE_data$logFC_abs <- abs(DE_data$logFC) # get absolute values
DE_data_abs <- subset(DE_data, select = c("gene_id", "logFC_abs"))


fwrite(gene_index, "BRCA/STN_filt/data/gene_index.tsv", col.names = F, sep = "\t")
fwrite(edge_list_index, "BRCA/STN_filt/data/edge_list_index.tsv", col.names = F, sep = "\t")
fwrite(DE_data_abs, "BRCA/STN_filt/data/updated_DE/logFC_scores_abs.tsv", col.names = F, sep = "\t")

