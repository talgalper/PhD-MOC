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

dir.create("results")

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# set seed for reproducability
set.seed(69)

# read in STRING data
weighted_string_data <- read.csv("data/weighted_string_data.csv", row.names = 1)
unweighted_string_data <- read.csv("data/unweighted_string_data.csv", row.names = 1)

# construct interactome
weighted_ppi <- construct_interactome(weighted_string_data)
unweighted_ppi <- construct_interactome(unweighted_string_data)

# read in seg mean data
weighted_data <- read.csv("data/weighted_pcsf_score.csv", row.names = 1)
unweighted_data <- read.csv("data/unweighted_pcsf_score.csv", row.names = 1)

# set terminals
weighted_terminals <- setNames(as.numeric(weighted_data$avg_seg_mean), weighted_data$ensembl_peptide_id)
unweighted_terminals <- setNames(as.numeric(unweighted_data$avg_seg_mean), unweighted_data$ensembl_peptide_id)

# run PCSF with random noise
weighted_subnet <- PCSF_rand(weighted_ppi, weighted_terminals, n = 50, r = 0.1, w = 2, b = 3, mu = 0.0005)
plot.PCSF(weighted_subnet, node_label_cex = 15)

unweighted_subnet <- PCSF_rand(unweighted_ppi, unweighted_terminals, n = 10, r = 0.1, w = 2, b = 3, mu = 0.0005)
plot.PCSF(unweighted_subnet, node_label_cex = 15)

# extract cluster data
weighted_clust <- clusters(weighted_subnet)
weighted_df <- data.frame(ensembl_peptide_id = names(weighted_clust$membership), cluster = factor(weighted_clust$membership))

# Calculate betweenness and centrality scores
weighted_betweenness <- betweenness(weighted_subnet) # number of shortest paths through node. number of times node acts as a bridge to other nodes
weighted_centrality <- degree(weighted_subnet) # how many connections each protein has

# Add betweenness and centrality scores to the cluster data frame
weighted_df$betweenness <- weighted_betweenness[as.character(weighted_df$ensembl_peptide_id)]
weighted_df$degree_centrality <- weighted_centrality[as.character(weighted_df$ensembl_peptide_id)]

# convert betweeness and centrality scores to integers
weighted_df$betweenness <- as.integer(weighted_df$betweenness)
weighted_df$degree_centrality <- as.integer(weighted_df$degree_centrality)

# format row names
rownames(weighted_df) <- 1:nrow(weighted_df)

# add a uniprot id column
weighted_uniprots <- getBM(attributes = c("ensembl_peptide_id", "uniprotswissprot"), 
                              filters = "ensembl_peptide_id", 
                              values = weighted_df$ensembl_peptide_id, 
                              mart = ensembl)

weighted_df <- merge(weighted_df, weighted_uniprots, by.x = "ensembl_peptide_id")

# format columns
weighted_df <- weighted_df[, c("ensembl_peptide_id", "uniprotswissprot", "cluster", "betweenness","degree_centrality")]
colnames(weighted_df)[2] <- "uniprot_id"

# sort by cluster number
weighted_df <- weighted_df[order(weighted_df$cluster), ]

# save to results dir
write.csv(weighted_df, "results/weighted_pcsf_results.csv", row.names = F)

# repeat for unweighted data
unweighted_clust <- clusters(unweighted_subnet)
unweighted_df <- data.frame(ensembl_peptide_id = names(unweighted_clust$membership), cluster = factor(unweighted_clust$membership))
unweighted_betweenness <- betweenness(unweighted_subnet) 
unweighted_centrality <- degree(unweighted_subnet) 
unweighted_df$betweenness <- unweighted_betweenness[as.character(unweighted_df$ensembl_peptide_id)]
unweighted_df$degree_centrality <- unweighted_centrality[as.character(unweighted_df$ensembl_peptide_id)]
unweighted_df$betweenness <- as.integer(unweighted_df$betweenness)
unweighted_df$degree_centrality <- as.integer(unweighted_df$degree_centrality)

rownames(unweighted_df) <- 1:nrow(unweighted_df)
unweighted_uniprots <- getBM(attributes = c("ensembl_peptide_id", "uniprotswissprot"), 
                           filters = "ensembl_peptide_id", 
                           values = unweighted_df$ensembl_peptide_id, 
                           mart = ensembl)
unweighted_df <- merge(unweighted_df, unweighted_uniprots, by.x = "ensembl_peptide_id")
unweighted_df <- unweighted_df[, c("ensembl_peptide_id", "uniprotswissprot", "cluster", "betweenness","degree_centrality")]
colnames(unweighted_df)[2] <- "uniprot_id"
unweighted_df <- unweighted_df[order(unweighted_df$cluster), ]

write.csv(unweighted_df, "results/unweighted_pcsf_results.csv", row.names = F)

