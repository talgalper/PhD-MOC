# check BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# required packages
required_packages <- c("biomaRt")

# Check if the required packages are installed, if not then install them
for (package in required_packages) {
  if (!requireNamespace(package, quietly = TRUE)) {
    BiocManager::install(package)
    library(package, character.only = TRUE)
  } else {
    library(package, character.only = TRUE)
  }
}

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")


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

weighted_df <- subset(merged_df, select = c("gene_id", "avg_seg_mean"))



unweighted_data <- read.csv("results/avg_unweighted_seg_mean.csv")
unweighted_data <- subset(unweighted_data, select = c("gene_id", "avg_seg_mean"))
colnames(unweighted_data) <- c("gene_id", "avg_seg_mean_unweighted")

unweighted_data <- merge(unweighted_data, weighted_df, by = "gene_id")
unweighted_data <- subset(unweighted_data, select = c("gene_id", "avg_seg_mean_unweighted"))
colnames(unweighted_data) <- c("gene_id", "avg_seg_mean")


write.csv(weighted_df, "results/avg_weighted_seg_mean_subset.csv")
write.csv(weighted_df, "results/avg_unweighted_seg_mean_subset.csv")




