
### requires fpocket_druggability.csv from druggability package run ###

weighted_data <- read.csv("results/avg_weighted_seg_mean.csv")
druggability_data <- read.csv("fpocket_druggability.csv")

gene_ids <- weighted_data$gene_id


# get the corresponding protein ensembles for each gene ID
uniprot_ids <- getBM(attributes = c("hgnc_symbol", "uniprotswissprot"), 
                     filters = "hgnc_symbol", 
                     values = gene_ids, 
                     mart = ensembl)

# merge the protein IDs with the original data frame based on the gene ID
new_df <- merge(weighted_data, uniprot_ids, by.x = "gene_id", by.y = "hgnc_symbol")
new_df <- subset(new_df, uniprotswissprot != "")

new_df <- new_df[, c("gene_id", "uniprotswissprot", "avg_seg_mean")]
colnames(new_df) <- c("gene_id", "uniprot_id", "avg_seg_mean")



merged_df <- merge(new_df, druggability_data, by = "uniprot_id")
merged_df <- merged_df[, c("gene_id", "uniprot_id", "avg_seg_mean", "pocket", "druggability", "num_drug_pockets")]
merged_df <- subset(merged_df, avg_seg_mean !="")

merged_df <- merged_df[merged_df$druggability >= 0.6, ]

final_df <- subset(merged_df, select = c("gene_id", "avg_seg_mean"))



