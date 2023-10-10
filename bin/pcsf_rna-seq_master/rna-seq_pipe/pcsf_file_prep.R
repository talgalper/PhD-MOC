library(tidyverse)
library(AnnotationDbi)
library(PCSF)
library(org.Hs.eg.db)
library(topGO)

tcga_dif_exp_data = "attemp_3_data/tcga_stage_II_vs_III_edgeR_dif_exp(p=0.1).tsv"
gene_list_output_filepath = "string_data/tcga_dif_exp_gene_list(p=0.1).tsv"
string_data_filepath = "attemp_3_data/STRING_network_edge.tsv"


# create gene list to download interaction data from STRING
data <- read.table(tcga_dif_exp_data, header = T, sep = "\t")
data$gene_id <- sub('\\.[0-9]*$', '', data$gene_id)

gene_list <- subset(data, select = c("gene_id"))

write_tsv(gene_list, gene_list_output_filepath, col_names = F)

# read in the data downloaded from STRING
string_data <- read.table(string_data_filepath, sep = "\t")
#string_data <- subset(string_data, select = c("V1", "V2", "V13"))


### initiate PCSF

# subset data and convert ensemble to gene symbol
terminals <- subset(data, select = c("gene_id", "logFC"))
terminals$symbol <- mapIds(org.Hs.eg.db, keys = terminals$gene_id, keytype = "ENSEMBL", column = "SYMBOL")
terminals <- subset(terminals, select = c("symbol", "logFC"))

# remove NA values if there are any
terminals <- na.omit(terminals)

# need to turn the list into a double
terminals <- setNames(c(terminals$logFC), c(terminals$symbol))

# run PCSF
ppi <- construct_interactome(string_data)
subnet <- PCSF(ppi, terminals, w = 2, b = 1, mu = 0.0005)

plot.PCSF(subnet, node_label_cex = 15)

# add random noise

subnet <- PCSF_rand(ppi, terminals, w = 2, b = 1, mu = 0.0005)

plot.PCSF(subnet)



