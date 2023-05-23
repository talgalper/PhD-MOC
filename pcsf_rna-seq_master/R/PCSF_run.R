devtools::install_github("IOR-Bioinformatics/PCSF", dependencies=TRUE, type="source", force=TRUE)

# check BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Check if the required packages are installed, if not then install them
if (!requireNamespace(package, quietly = TRUE)) {
    BiocManager::install("biomaRt")
}

library(PCSF)
library(biomaRt)


ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# empty global environment except for ensembl object
rm(list = ls()[!ls() %in% c("ensembl")])


# set seed for reproducability
set.seed(69)

stage <- "stage_I"

# read in STRING data
tcga_string_data <- read.csv(paste0("pcsf/", stage, "/", "tcga/tcga_string_data.csv"), row.names = 1)
kylie_string_data <- read.csv(paste0("pcsf/", stage, "/", "kylie/kylie_string_data.csv"), row.names = 1)

# construct interactome
tcga_ppi <- construct_interactome(tcga_string_data)
kylie_ppi <- construct_interactome(kylie_string_data)

# read in score file
tcga_data <- read.csv(paste0("pcsf/", stage, "/", "tcga/pcsf_tcga_score.csv"), row.names = 1)
kylie_data <- read.csv(paste0("pcsf/", stage, "/", "kylie/pcsf_kylie_score.csv"), row.names = 1)

# set terminals
tcga_terminals <- setNames(as.numeric(tcga_data$logFC), tcga_data$uniprotswissprot)
kylie_terminals <- setNames(as.numeric(kylie_data$logFC), kylie_data$uniprotswissprot)

# run PCSF with random noise
tcga_subnet <- PCSF_rand(tcga_ppi, tcga_terminals, n = 50, r = 0.1, w = 2, b = 3, mu = 0.0005)
plot.PCSF(tcga_subnet, node_label_cex = 15)

kylie_subnet <- PCSF_rand(kylie_ppi, kylie_terminals, n = 10, r = 0.1, w = 2, b = 3, mu = 0.0005)
plot.PCSF(kylie_subnet, node_label_cex = 15)

# extract cluster data
tcga_clust <- clusters(tcga_subnet)
tcga_df <- data.frame(uniprot_id = names(tcga_clust$membership), cluster = factor(tcga_clust$membership))

# Calculate betweenness and centrality scores
tcga_betweenness <- betweenness(tcga_subnet) # number of shortest paths through node. number of times node acts as a bridge to other nodes
tcga_centrality <- degree(tcga_subnet) # how many connections each protein has

# Add betweenness and centrality scores to the cluster data frame
tcga_df$betweenness <- tcga_betweenness[as.character(tcga_df$uniprot_id)]
tcga_df$degree_centrality <- tcga_centrality[as.character(tcga_df$uniprot_id)]

# convert betweeness and centrality scores to integers
tcga_df$betweenness <- as.integer(tcga_df$betweenness)
tcga_df$degree_centrality <- as.integer(tcga_df$degree_centrality)

# format row names
rownames(tcga_df) <- 1:nrow(tcga_df)

## add a uniprot id column
#weighted_uniprots <- getBM(attributes = c("ensembl_peptide_id", "uniprotswissprot"), 
#                              filters = "ensembl_peptide_id", 
#                              values = weighted_df$ensembl_peptide_id, 
#                              mart = ensembl)
#
#weighted_df <- merge(weighted_df, weighted_uniprots, by.x = "ensembl_peptide_id")
#
## format columns
#weighted_df <- weighted_df[, c("ensembl_peptide_id", "uniprotswissprot", "cluster", "betweenness","degree_centrality")]
#colnames(weighted_df)[2] <- "uniprot_id"

# sort by cluster number
tcga_df <- tcga_df[order(tcga_df$cluster), ]

# save to results dir
write.csv(tcga_df, paste0("pcsf/", stage, "/", "tcga/results/tcga_pcsf_results.csv"), row.names = F)

# repeat for kylie data
kylie_clust <- clusters(kylie_subnet)
kylie_df <- data.frame(uniprot_id = names(kylie_clust$membership), cluster = factor(kylie_clust$membership))
kylie_betweenness <- betweenness(kylie_subnet) 
kylie_centrality <- degree(kylie_subnet) 
kylie_df$betweenness <- kylie_betweenness[as.character(kylie_df$uniprot_id)]
kylie_df$degree_centrality <- kylie_centrality[as.character(kylie_df$uniprot_id)]
kylie_df$betweenness <- as.integer(kylie_df$betweenness)
kylie_df$degree_centrality <- as.integer(kylie_df$degree_centrality)

rownames(kylie_df) <- 1:nrow(kylie_df)

#unweighted_uniprots <- getBM(attributes = c("ensembl_peptide_id", "uniprotswissprot"), 
#                           filters = "ensembl_peptide_id", 
#                           values = unweighted_df$ensembl_peptide_id, 
#                           mart = ensembl)
#unweighted_df <- merge(unweighted_df, unweighted_uniprots, by.x = "ensembl_peptide_id")
#unweighted_df <- unweighted_df[, c("ensembl_peptide_id", "uniprotswissprot", "cluster", "betweenness","degree_centrality")]
#colnames(unweighted_df)[2] <- "uniprot_id"

kylie_df <- kylie_df[order(kylie_df$cluster), ]

write.csv(kylie_df, paste0("pcsf/", stage, "/", "kylie/results/kylie_pcsf_results.csv"), row.names = F)

