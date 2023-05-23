
# check BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# required packages
required_packages <- c("biomaRt", "dplyr", "tidyr", "RCy3")

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


# empty global environment except for ensembl object
rm(list = ls()[!ls() %in% c("ensembl")])

# identify which stage you want to run
stage <- "stage_IV"

# read in the data
tcga_data <- read.csv(paste0("edgeR_results/tcga_data/", stage, "/", stage, "_edgeR_hits.csv"))
kylie_data <- read.csv(paste0("edgeR_results/kylie_data/", stage, "/", stage, "_edgeR_hits.csv"))

# subset the gene_id and logFC columns
tcga_data <- subset(tcga_data, select = c("gene_id", "logFC"))
kylie_data <- subset(kylie_data, select = c("gene_id", "logFC"))

hist(kylie_data$logFC)
hist(tcga_data$logFC)

# remove all values with logFC between -1 and 1
tcga_data$logFC <- ifelse(tcga_data$logFC > -1 & tcga_data$logFC < 1, NA, tcga_data$logFC)
tcga_data <- na.omit(tcga_data)

kylie_data$logFC <- ifelse(kylie_data$logFC > -1 & kylie_data$logFC < 1, NA, kylie_data$logFC)
kylie_data <- na.omit(kylie_data)

# convert to protein ensembl
tcga_uniprot <- getBM(attributes = c("ensembl_gene_id", "uniprotswissprot"), 
                      filters = "ensembl_gene_id", 
                      values = tcga_data$gene_id, 
                      mart = ensembl)

kylie_uniprot <- getBM(attributes = c("ensembl_gene_id", "uniprotswissprot"), 
                       filters = "ensembl_gene_id", 
                       values = kylie_data$gene_id, 
                       mart = ensembl)


# remove empty rows
tcga_uniprot <- subset(tcga_uniprot, uniprotswissprot != "")
kylie_uniprot <- subset(kylie_uniprot, uniprotswissprot != "")

# check for duplicate protein ensembles
tcga_uniprot <- distinct(tcga_uniprot)
kylie_uniprot <- distinct(kylie_uniprot)

# merge back with original data
colnames(tcga_uniprot)[1] <- "gene_id"
colnames(kylie_uniprot)[1] <- "gene_id"

tcga_protein_data <- merge(tcga_uniprot, tcga_data, by = "gene_id")
kylie_protein_data <- merge(kylie_uniprot, kylie_data, by = "gene_id")

# create score files
tcga_protein_data <- subset(tcga_protein_data, select = c("uniprotswissprot", "logFC"))
kylie_protein_data <- subset(kylie_protein_data, select = c("uniprotswissprot", "logFC"))

write.csv(tcga_protein_data, paste0("pcsf/", stage, "/", "tcga/pcsf_tcga_score.csv"))
write.csv(kylie_protein_data, paste0("pcsf/", stage, "/", "kylie/pcsf_kylie_score.csv"))


## create ensemble protein list for cytoscape
write.table(tcga_protein_data$uniprotswissprot, paste0("pcsf/", stage, "/", "tcga/cytoscape/tcga_proteins.txt"), row.names = F, col.names = F, quote = F)
write.table(kylie_protein_data$uniprotswissprot, paste0("pcsf/", stage, "/", "kylie/cytoscape/kylie_proteins.txt"), row.names = F, col.names = F, quote = F)

## format cytoscape output
string_edge_data <- read.table(paste0("pcsf/", stage, "/", "tcga/cytoscape/STRING tcga network default edge.csv"), header = T, sep = ",", stringsAsFactors = F)
ppi_list <- subset(string_edge_data, select = c("name", "stringdb..score"))
ppi_list <- ppi_list %>% 
  separate(name, sep = " ", into = c("node_1", "del", "node_2"))
ppi_list <- subset(ppi_list, select = c("node_1", "node_2", "stringdb..score"))
ppi_list$node_1 <- gsub(".*.\\.", "", ppi_list$node_1)
ppi_list$node_2 <- gsub(".*.\\.", "", ppi_list$node_2)

## replace protein ensembles with uniprot IDs
string_node_data <- read.table(paste0("pcsf/", stage, "/", "tcga/cytoscape/STRING tcga network default node.csv"), header = T, sep = ",", stringsAsFactors = F)
node_list <- subset(string_node_data, select = c("name", "query.term"))
node_list$name <- gsub(".*.\\.", "", node_list$name)
# Add an additional column to preserve the original row order
ppi_list$original_order <- seq_len(nrow(ppi_list))
# Merge the ppi_list and node_list
merged_df <- merge(ppi_list, node_list, by.x = "node_1", by.y = "name", all.x = TRUE)
merged_df <- merge(merged_df, node_list, by.x = "node_2", by.y = "name", all.x = TRUE)
# Sort the merged data frame based on the original order 
merged_df <- merged_df[order(merged_df$original_order), ]
# Select the required columns
final_df <- merged_df[, c("query.term.x", "query.term.y", "stringdb..score")]
# Rename the columns
colnames(final_df) <- c("node_1", "node_2", "uniprot_id")


write.csv(final_df, paste0("pcsf/", stage, "/", "tcga/tcga_string_data.csv"))

# repeat for kylie data
string_edge_data <- read.table(paste0("pcsf/", stage, "/", "kylie/cytoscape/STRING kylie network default edge.csv"), header = T, sep = ",", stringsAsFactors = F)
ppi_list <- subset(string_edge_data, select = c("name", "stringdb..score"))
ppi_list <- ppi_list %>% 
  separate(name, sep = " ", into = c("node_1", "del", "node_2"))
ppi_list <- subset(ppi_list, select = c("node_1", "node_2", "stringdb..score"))
ppi_list$node_1 <- gsub(".*.\\.", "", ppi_list$node_1)
ppi_list$node_2 <- gsub(".*.\\.", "", ppi_list$node_2)

# replace protein ensembles with uniprot IDs
string_node_data <- read.table(paste0("pcsf/", stage, "/", "kylie/cytoscape/STRING kylie network default node.csv"), header = T, sep = ",", stringsAsFactors = F)
node_list <- subset(string_node_data, select = c("name", "query.term"))
node_list$name <- gsub(".*.\\.", "", node_list$name)
ppi_list$original_order <- seq_len(nrow(ppi_list))
merged_df <- merge(ppi_list, node_list, by.x = "node_1", by.y = "name", all.x = TRUE)
merged_df <- merge(merged_df, node_list, by.x = "node_2", by.y = "name", all.x = TRUE)
merged_df <- merged_df[order(merged_df$original_order), ]
final_df <- merged_df[, c("query.term.x", "query.term.y", "stringdb..score")]
colnames(final_df) <- c("node_1", "node_2", "uniprot_id")

write.csv(final_df, paste0("pcsf/", stage, "/", "kylie/kylie_string_data.csv"))

