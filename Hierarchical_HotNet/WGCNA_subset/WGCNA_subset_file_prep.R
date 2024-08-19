library(tidyverse)

# load in BRCA data
load("../BRCA_pipe/latest_run/RData/STRING_PPI_FULL.RData")
load("../WGCNA/BRCA/RData/DE_subset/dif_exp.RData")
DE_data <- subset(dif_exp, select = c("gene_id", "logFC"))
DE_data$logFC_abs <- abs(DE_data$logFC) # get absolute values

load("../WGCNA/BRCA/RData/all_default/signed/venn_data.RData")
common_genes <- Reduce(intersect, list(DE_genes, tumour_associated, top_kwithin, top_gene_membership))
DE_data <- DE_data[DE_data$gene_id %in% common_genes, ]

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

# remove interaction scores
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


write_tsv(gene_index, "WGCNA_subset/data/gene_index.tsv", col_names = F)
write_tsv(edge_list_index, "WGCNA_subset/data/edge_list_index.tsv", col_names = F)
write_tsv(final_df, "WGCNA_subset/data/STRING_network.tsv", col_names = F)
write_tsv(DE_data, "WGCNA_subset/data/logFC_scores.tsv", col_names = F)
write_tsv(DE_data_abs, "WGCNA_subset/data/logFC_scores_abs.tsv", col_names = F)

