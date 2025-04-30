library(data.table)

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
DE_data <- read.csv("../MOC_pipe/DE/MOC_vs_BEN/MOC_DE_results.csv")
DE_data <- subset(DE_data, select = c("ensembl_gene_id", "logFC"))
DE_data$logFC_abs <- abs(DE_data$logFC) # get absolute values
DE_data_abs <- subset(DE_data, select = c("ensembl_gene_id", "logFC_abs"))


fwrite(gene_index, "MOC/data/gene_index.tsv", col.names = F, sep = "\t")
fwrite(edge_list_index, "MOC/data/edge_list_index.tsv", col.names = F, sep = "\t")
fwrite(DE_data_abs, "MOC/data/logFC_scores_abs.tsv", col.names = F, sep = "\t")




# create score files for MOC vs GTEx
DE_data <- read.csv("../MOC_pipe/DE/MOC_vs_GTEx/MOC_vs_GTEx_DE_results.csv")
DE_data <- subset(DE_data, select = c("ensembl_gene_id", "logFC"))
DE_data$logFC_abs <- abs(DE_data$logFC) # get absolute values
DE_data_abs <- subset(DE_data, select = c("ensembl_gene_id", "logFC_abs"))

fwrite(gene_index, "MOC/MOC_vs_GTEx/data/gene_index.tsv", col.names = F, sep = "\t")
fwrite(edge_list_index, "MOC/MOC_vs_GTEx/data/edge_list_index.tsv", col.names = F, sep = "\t")
fwrite(DE_data_abs, "MOC/MOC_vs_GTEx/data/logFC_scores_abs.tsv", col.names = F, sep = "\t")







