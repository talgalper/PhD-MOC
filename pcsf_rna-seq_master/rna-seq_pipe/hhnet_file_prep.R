### generate all the files necessary to run Hierarchical Hotnet


library(tidyverse)
library(org.Hs.eg.db)


input_file = "stage_I_master_df.tsv"
stage <- "stage_I"
output_dir <- "hhnet"




## only required to be run once
# create list of genes for STRING to collect interaction scores
gene_list <- read.table(input_file, header = T, sep = "\t")
gene_list <- as.data.frame(index_to_gene$gene_id)
write_tsv(gene_list, "hhnet/string/gene_list.tsv", col_names = F)

# read in STRING interaction data and format
string_data <- read.table("hhnet/string/STRING_network_edge.csv", header = T, sep = ",")
ppi_list <- subset(string_data, select = c("name"))
ppi_list <- ppi_list %>% 
  separate(name, sep = " ", into = c("node_1", "del", "node_2"))
ppi_list <- subset(ppi_list, select = c("node_1", "node_2"))
ppi_list$node_1 <- gsub(".*.\\.", "", ppi_list$node_1)
ppi_list$node_2 <- gsub(".*.\\.", "", ppi_list$node_2)





# gene to score file (tcga data)
gene_to_score <- read.table(input_file, header = T, sep = "\t")
gene_to_score$gene_id <- sub('\\.[0-9]*$', '', gene_to_score$gene_id)
#gene_to_score <- subset(gene_to_score, select = c("gene_id", "logFC"))
#gene_to_score$gene_id <- mapIds(org.Hs.eg.db, keys = gene_to_score$gene_id, keytype = "ENSEMBL", column = "SYMBOL")

# remove NA values if there are any
gene_to_score <- na.omit(gene_to_score)

write_tsv(gene_to_score, file = file.path(output_dir, stage, "gene_to_score.tsv"), col_names = F)


# index to gene
index_to_gene <- subset(gene_to_score, select = c("gene_id"))
index_to_gene$row_num <- seq.int(nrow(index_to_gene))
index_to_gene <- index_to_gene[c("row_num", "gene_id")]

write_tsv(index_to_gene, file = file.path(output_dir, stage, "index_to_gene.tsv"), col_names = F)



# create edge list for hhnet (indexed edge list). This code is from Sam's github.
edge_list <- as.matrix(ppi_list)

graph_index <- edge_list %>%
  as.vector() %>%
  unique() %>%
  sort() %>%
  data.frame(name = ., index = 1:length(.))

graph_index_vec <- setNames(
  graph_index$index, 
  graph_index$name
)

edge_list_index <- matrix(graph_index_vec[edge_list], ncol = 2)
edge_list_index <- as.data.frame(edge_list_index)


write_tsv(edge_list_index, "hhnet/edge_list.tsv")




