
# check BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# required packages
required_packages <- c("biomaRt", "dplyr", "tidyr")

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

### create PCSF files ###

weighted_means <- read.csv("data/cnv_pipe/avg_weighted_seg_mean.csv", row.names = 1)
unweighted_means <- read.csv("data/cnv_pipe/avg_unweighted_seg_mean.csv", row.names = 1)

hist(weighted_means$avg_seg_mean, xlim = c(1, -1))
hist(unweighted_means$avg_seg_mean, xlim = c(1, -1))

# remove all values with seg_mean between -0.3 and 0.3
weighted_means$avg_seg_mean <- ifelse(weighted_means$avg_seg_mean > -0.3 & weighted_means$avg_seg_mean < 0.3, NA, weighted_means$avg_seg_mean)
weighted_means <- na.omit(weighted_means)

unweighted_means$avg_seg_mean <- ifelse(unweighted_means$avg_seg_mean > -0.3 & unweighted_means$avg_seg_mean < 0.3, NA, unweighted_means$avg_seg_mean)
unweighted_means <- na.omit(unweighted_means)

weighted_gene_ids <- weighted_means$gene_id
unweighted_gene_ids <- unweighted_means$gene_id

# convert to protein ensembl
weighted_protein_ens <- getBM(attributes = c("hgnc_symbol", "ensembl_peptide_id"), 
                     filters = "hgnc_symbol", 
                     values = weighted_gene_ids, 
                     mart = ensembl)

unweighted_protein_ens <- getBM(attributes = c("hgnc_symbol", "ensembl_peptide_id"), 
                              filters = "hgnc_symbol", 
                              values = unweighted_gene_ids, 
                              mart = ensembl)

# remove empty rows
weighted_ens <- subset(weighted_protein_ens, ensembl_peptide_id != "")
unweighted_ens <- subset(unweighted_protein_ens, ensembl_peptide_id != "")

# check for duplicate protein ensembles
weighted_ens <- distinct(weighted_ens)
unweighted_ens <- distinct(unweighted_ens)

# merge back with original data
weighted_ens_means <- merge(weighted_means, weighted_ens, by.x = "gene_id", by.y = "hgnc_symbol")
unweighted_ens_means <- merge(unweighted_means, unweighted_ens, by.x = "gene_id", by.y = "hgnc_symbol")

# reorder columns
weighted_ens_means <- weighted_ens_means[, c("gene_id", "ensembl_peptide_id", "avg_seg_mean")]
unweighted_ens_means <- unweighted_ens_means[, c("gene_id", "ensembl_peptide_id", "avg_seg_mean")]

# create score files
weighted_data <- subset(weighted_ens_means, select = c("ensembl_peptide_id", "avg_seg_mean"))
unweighted_data <- subset(unweighted_ens_means, select = c("ensembl_peptide_id", "avg_seg_mean"))

write.csv(weighted_data, "data/weighted_pcsf_score.csv")
write.csv(unweighted_data, "data/unweighted_pcsf_score.csv")


## create ensemble protein list for cytoscape
write.table(weighted_ens_means$ensembl_peptide_id, "data/cytoscape_stuff/weighted_proteins.txt", row.names = F, col.names = F, quote = F)
write.table(unweighted_ens_means$ensembl_peptide_id, "data/cytoscape_stuff/unweighted_proteins.txt", row.names = F, col.names = F, quote = F)

## format cytoscape output
# weighted PPI list
string_data <- read.table("data/cytoscape_stuff/STRING network weighted edge.csv", header = T, sep = ",", stringsAsFactors = F)
ppi_list <- subset(string_data, select = c("name", "stringdb..score"))
ppi_list <- ppi_list %>% 
  separate(name, sep = " ", into = c("node_1", "del", "node_2"))
ppi_list <- subset(ppi_list, select = c("node_1", "node_2", "stringdb..score"))
ppi_list$node_1 <- gsub(".*.\\.", "", ppi_list$node_1)
ppi_list$node_2 <- gsub(".*.\\.", "", ppi_list$node_2)

write.csv(ppi_list, "data/weighted_string_data.csv")

# unweighted PPI list
string_data <- read.table("data/cytoscape_stuff/STRING network unweighted edge.csv", header = T, sep = ",", stringsAsFactors = F)
ppi_list <- subset(string_data, select = c("name", "stringdb..score"))
ppi_list <- ppi_list %>% 
  separate(name, sep = " ", into = c("node_1", "del", "node_2"))
ppi_list <- subset(ppi_list, select = c("node_1", "node_2", "stringdb..score"))
ppi_list$node_1 <- gsub(".*.\\.", "", ppi_list$node_1)
ppi_list$node_2 <- gsub(".*.\\.", "", ppi_list$node_2)

write.csv(ppi_list, "data/unweighted_string_data.csv")
