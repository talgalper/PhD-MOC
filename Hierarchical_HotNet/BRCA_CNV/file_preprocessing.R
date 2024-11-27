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
load("../BRCA_CNV/RData/BRCA_CNV_alterationsPerc.RData")
cnv_summary_filt <- subset(cnv_summary_filt, select = c("ensembl_gene_id", "alterations_percentage"))
cnv_summary_filt <- cnv_summary_filt[cnv_summary_filt$ensembl_gene_id %in% gene_index$ensembl_id, ] # keep genes present in both PPI and CNV
cnv_summary_filt$alterations_percentage <- cnv_summary_filt$alterations_percentage / 100 # convert to 0-1

fwrite(gene_index, "BRCA_CNV/data/gene_index.tsv", col.names = F, sep = "\t")
fwrite(edge_list_index, "BRCA_CNV/data/edge_list_index.tsv", col.names = F, sep = "\t")
fwrite(cnv_summary_filt, "BRCA_CNV/data/cnv_alterations.tsv", col.names = F, sep = "\t")

