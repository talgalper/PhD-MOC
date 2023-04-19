library(biomaRt)
library(tidyr)
library(readr)


### create protein (ensembl) to score file

# load ensembl
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# load weighted copy number data
weighted_data <- read.csv("weighted_seg_mean_466074.csv", row.names = 1)
df <- weighted_data # save to new variable as backup

# create a vector of gene IDs from your data frame
gene_ids <- unique(df$gene_id)

# get the corresponding protein ensembles for each gene ID
protein_ens <- getBM(attributes = c("hgnc_symbol", "ensembl_peptide_id"), 
                     filters = "hgnc_symbol", 
                     values = gene_ids, 
                     mart = ensembl)

# merge the protein IDs with the original data frame based on the gene ID
protein_to_score <- merge(df, protein_ens, by.x = "gene_id", by.y = "hgnc_symbol")

# drop the gene_id column and remove empty rows
protein_to_score <- subset(protein_to_score, ensembl_peptide_id != "")
protein_to_score$gene_id <- NULL

protein_to_score <- protein_to_score[, c("ensembl_peptide_id", "seg_mean")]

# seperate + and - values
positive_df <- protein_to_score[protein_to_score$seg_mean > 0, ]
negative_df <- protein_to_score[protein_to_score$seg_mean <= 0, ]

# convert all values in negative_df to positive
negative_df$seg_mean <- abs(negative_df$seg_mean)


### create a indexed edge list file for + and - ###

## cytoscape input
# create list of proteins for cytoscape
protein_list <- protein_to_score$ensembl_peptide_id
protein_list <- as.data.frame(protein_list)
write.table(gene_list, "protein_list.txt", col.names = F, row.names = F, quote = F)


## cytoscape output
string_data <- read.table("single_run_466074/STRING protein ensembl default edge.csv", header = T, sep = ",", stringsAsFactors = F)
ppi_list <- subset(string_data, select = c("name"))
ppi_list <- ppi_list %>% 
  separate(name, sep = " ", into = c("node_1", "del", "node_2"))
ppi_list <- subset(ppi_list, select = c("node_1", "node_2"))
ppi_list$node_1 <- gsub(".*.\\.", "", ppi_list$node_1)
ppi_list$node_2 <- gsub(".*.\\.", "", ppi_list$node_2)


## Compare proteins in positive_df and ppi_list, and keep only those identified in both
upreg_common_proteins <- intersect(positive_df$ensembl_peptide_id, ppi_list$node_1)
positive_df <- positive_df[positive_df$ensembl_peptide_id %in% upreg_common_proteins, ]
upreg_ppi_list <- ppi_list[ppi_list$node_1 %in% upreg_common_proteins, ]

# Repeat the above for negative_df
dnreg_common_proteins <- intersect(negative_df$ensembl_peptide_id, ppi_list$node_1)
negative_df <- negative_df[negative_df$ensembl_peptide_id %in% dnreg_common_proteins, ]
dnreg_ppi_list <- ppi_list[ppi_list$node_1 %in% dnreg_common_proteins, ]


write_tsv(positive_df, "hierarchical-hotnet_upreg/data/upreg_protein_to_score.tsv", col_names = F)
write_tsv(negative_df, "hierarchical-hotnet_dnreg/data/dnreg_protein_to_score.tsv", col_names = F)


### create protein index file ###
upreg_index_to_protein <- subset(positive_df, select = c("ensembl_peptide_id"))
upreg_index_to_protein$row_num <- seq.int(nrow(upreg_index_to_protein))
upreg_index_to_protein <- upreg_index_to_protein[c("row_num", "ensembl_peptide_id")]

dnreg_index_to_protein <- subset(negative_df, select = c("ensembl_peptide_id"))
dnreg_index_to_protein$row_num <- seq.int(nrow(dnreg_index_to_protein))
dnreg_index_to_protein <- dnreg_index_to_protein[c("row_num", "ensembl_peptide_id")]

write_tsv(upreg_index_to_protein, "hierarchical-hotnet_upreg/data/upreg_index_to_protein.tsv", col_names = F)
write_tsv(dnreg_index_to_protein, "hierarchical-hotnet_dnreg/data/dnreg_index_to_protein.tsv", col_names = F)


### create a indexed edge list file ###
upreg_edge_list_index<- data.frame(from = match(upreg_ppi_list$node_1, upreg_index_to_protein$ensembl_peptide_id),
                              to = match(upreg_ppi_list$node_2, upreg_index_to_protein$ensembl_peptide_id))
upreg_edge_list_index <- na.omit(upreg_edge_list_index)
write_tsv(upreg_edge_list_index, "hierarchical-hotnet_upreg/data/upreg_edge_list_index.tsv", col_names = F)


# create indexed edge list for downregulated proteins
dnreg_edge_list_index <- data.frame(from = match(dnreg_ppi_list$node_1, dnreg_index_to_protein$ensembl_peptide_id),
                              to = match(dnreg_ppi_list$node_2, dnreg_index_to_protein$ensembl_peptide_id))
dnreg_edge_list_index <- na.omit(dnreg_edge_list_index)
write_tsv(dnreg_edge_list_index, "hierarchical-hotnet_dnreg/data/dnreg_edge_list_index.tsv", col_names = F)



