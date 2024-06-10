library(TCGAbiolinks)
library(SummarizedExperiment)
library(edgeR)
library(tidyverse)
library(biomaRt)
library(PCSF)
library(progress)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")


load("RData/DE_results_master_paired.RData")
dif_exp <- DE_results$TCGA_lumB$dif_exp


data <- subset(dif_exp, select = c("gene_id", "logFC"))

# convert to gene symbol
gene_id <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), 
                 filters = "ensembl_gene_id", 
                 values = data$gene_id, 
                 mart = ensembl)

# remove empty rows
gene_id <- subset(gene_id, external_gene_name != "")

# check for duplicate uniprot ids
gene_id <- distinct(gene_id)

# merge back with original data
colnames(gene_id)[1] <- "gene_id"

gene_data <- merge(gene_id, data, by = "gene_id")

gene_data <- subset(gene_data, select = c("external_gene_name", "logFC"))

# check to see if genes that at all genes were converted at least once
missing_genes <- anti_join(data, gene_id, by = "gene_id")

write.table(gene_data$external_gene_name, "intermediate/paired/lumB/gene_list.txt", quote = F, row.names = F, col.names = F)


## create table with details of lost genes
missing_genes_convert <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "uniprot_gn_id", "description"), 
                               filters = "ensembl_gene_id", 
                               values = missing_genes$gene_id, 
                               mart = ensembl)
missing_genes_convert <- merge(missing_genes_convert, missing_genes, by.x = "ensembl_gene_id", by.y = "gene_id", all = T)

# number of unrecognised terms
table(is.na(missing_genes_convert$uniprot_gn_id) & is.na(missing_genes_convert$description))

novel_transcripts <- missing_genes_convert[grep("novel transcript", missing_genes_convert$description), ]
novel_proteins <- missing_genes_convert[grep("novel protein", missing_genes_convert$description), ]
pseudogene <- missing_genes_convert[grep("pseudogene", missing_genes_convert$description), ]



#### get interaction data ####
string_edge_data <- read.table("intermediate/paired/lumB/STRING network (physical) default edge.csv", header = T, sep = ",", stringsAsFactors = F)
ppi_list <- subset(string_edge_data, select = c("name", "stringdb..score"))
ppi_list <- ppi_list %>% 
  separate(name, sep = " ", into = c("node_1", "del", "node_2"))
ppi_list <- subset(ppi_list, select = c("node_1", "node_2", "stringdb..score"))
ppi_list$node_1 <- gsub(".*.\\.", "", ppi_list$node_1)
ppi_list$node_2 <- gsub(".*.\\.", "", ppi_list$node_2)

string_node_data <- read.table("intermediate/paired/lumB/STRING network (physical) default node.csv", header = T, sep = ",", stringsAsFactors = F)
node_list <- subset(string_node_data, select = c("name", "query.term"))
node_list$name <- gsub(".*.\\.", "", node_list$name)
ppi_list$original_order <- seq_len(nrow(ppi_list))
merged_df <- merge(ppi_list, node_list, by.x = "node_1", by.y = "name", all.x = TRUE)
merged_df <- merge(merged_df, node_list, by.x = "node_2", by.y = "name", all.x = TRUE)
merged_df <- merged_df[order(merged_df$original_order), ]

final_df <- merged_df[, c("query.term.x", "query.term.y", "stringdb..score")]
colnames(final_df) <- c("node_1", "node_2", "score")

save(final_df, gene_data, file = "RData/paired/lumB/PCSF_input.RData")


# set seed for reproducibility 
set.seed(1234)

# construct interactome
ppi <- construct_interactome(final_df)

# set terminals
terminals <- setNames(as.numeric(gene_data$logFC), gene_data$external_gene_name)

# run PCSF with random noise

# time a pcsf run
start_time <- Sys.time()
subnet <- PCSF_rand(ppi, terminals, n = 50, r = 0.1, w = 2, b = 1, mu = 0.0005)
elapsed_time <- Sys.time() - start_time
print(elapsed_time)

plot.PCSF(subnet, node_label_cex = 15)

save(subnet, file = "RData/paired/lumB/PCSF_subnet.RData")




# extract cluster data
clust <- components(subnet)
df <- data.frame(gene_id = names(clust$membership), cluster = factor(clust$membership))
betweenness <- betweenness(subnet) 
centrality <- degree(subnet) 
df$betweenness <- betweenness[as.character(df$gene_id)]
df$degree_centrality <- centrality[as.character(df$gene_id)]
df$betweenness <- as.integer(df$betweenness)
df$degree_centrality <- as.integer(df$degree_centrality)

rownames(df) <- 1:nrow(df)

df <- df[order(-df$degree_centrality), ]

write.csv(df, "intermediate/paired/lumB/PCSF_output.csv", row.names = F)

# convert external gene name to uniprot
gn_to_uniprot <- getBM(attributes = c("external_gene_name", "uniprot_gn_id"), 
                       filters = "external_gene_name", 
                       values = df$gene_id, 
                       mart = ensembl)

# get list of missing genes
missing_genes <- gn_to_uniprot[gn_to_uniprot$uniprot_gn_id == "", ]
missing_genes <- missing_genes$external_gene_name

PCSF_master <- merge(gn_to_uniprot, df, by.x = "external_gene_name", by.y = "gene_id")

PCSF_master <- merge(PCSF_master, gene_data, by = "external_gene_name")
PCSF_master <- PCSF_master[PCSF_master$uniprot_gn_id != "", ]


# load Fpocket data
af_drugability <- read.csv("../druggability_results/fpocket_druggability.csv")

# merge PCSF data with AF
PCSF_results <- merge(PCSF_master, af_drugability, by.x = "uniprot_gn_id", by.y = "uniprot_id")

missing_genes <- unique(PCSF_master$external_gene_name)[!unique(PCSF_master$external_gene_name) %in% unique(PCSF_results$external_gene_name)]

write.csv(PCSF_results, "intermediate/paired/LumB/PCSF_druggability.csv")



# keep structure duplicates with highest druggability
filtered_results <- PCSF_results %>%
  group_by(external_gene_name) %>%
  filter(druggability == max(druggability))

write.csv(filtered_results, "intermediate/paired/LumB/PCSF_master_unique.csv")


## cryptic pocket scores
pcsf_master <- filtered_results

pocketminer_data <- read.csv("../pocketminer/results/pocketminer_results_3.0.csv")
pcsf_master <- merge(pcsf_master, pocketminer_data, by = "ID")

missing_genes <- filtered_results$external_gene_name[!filtered_results$external_gene_name %in% pcsf_master$external_gene_name]


#### Ranking ####
betweeness_norm <- (pcsf_master$betweenness - min(pcsf_master$betweenness)) / (max(pcsf_master$betweenness) - min(pcsf_master$betweenness))

centrality_norm <- (pcsf_master$degree_centrality - min(pcsf_master$degree_centrality)) / (max(pcsf_master$degree_centrality) - min(pcsf_master$degree_centrality))

# Define the number of weight values (genes in top x) to test (including 0)
# i.e. 0, 0.1, 0.2, 0.3, ..., 1
num_weights <- 11

# Create weight value sequences that sum up to 1
betweeness_w_values <- seq(0, 1, length.out = num_weights)
centrality_w_values <- seq(0, 1, length.out = num_weights)
druggability_w_values <- seq(0, 1, length.out = num_weights)
cryptic_pocket_w_values <- seq(0, 1, length.out = num_weights)


# Calculate the total number of iterations a.k.a number of variables
total_iterations <- num_weights ^ 4

# Create a progress bar
pb <- progress_bar$new(total = total_iterations)

# Create an empty data frame with proper column names
sensitivity_results <- data.frame(
  betweeness_weight = numeric(),
  centrality_weight = numeric(),
  druggability_weight = numeric(),
  cryptic_pocket_weight = numeric(),
  top_genes = character(),
  stringsAsFactors = FALSE
)

# Perform sensitivity test
for (betweeness_w in betweeness_w_values) {
  for (centrality_w in centrality_w_values) {
    for (druggability_w in druggability_w_values) {
      for (cryptic_pocket_w in cryptic_pocket_w_values) {
        # Check if the weights sum up to 1
        if (betweeness_w + centrality_w + druggability_w + cryptic_pocket_w == 1) {
          # Combine scores using current weights
          combined_score <- (betweeness_w * betweeness_norm) +
            (centrality_w * centrality_norm) +
            (druggability_w * pcsf_master$druggability) +
            (cryptic_pocket_w * pcsf_master$max_hit)
          
          # Rank the genes based on the combined score
          pcsf_master_ranked <- pcsf_master
          pcsf_master_ranked$combined_score <- combined_score
          pcsf_master_ranked <- pcsf_master_ranked[order(-pcsf_master_ranked$combined_score), ]
          
          # Select the top 10 genes
          top_genes <- head(pcsf_master_ranked$external_gene_name, 10)
          
          # Add results to the sensitivity_results data frame
          sensitivity_results <- rbind(sensitivity_results, list(
            betweeness_weight = betweeness_w,
            centrality_weight = centrality_w,
            druggability_weight = druggability_w,
            cryptic_pocket_weight = cryptic_pocket_w,
            top_genes = paste(top_genes, collapse = ', ')
          ))
        }
        
        # Increment the progress bar
        pb$tick()
      }
    }
  }
}


# Split the "top_genes" column into a list of genes
sensitivity_results$top_genes_list <- strsplit(sensitivity_results$top_genes, ', ')

# Count the occurrences of each gene
all_genes <- unlist(sensitivity_results$top_genes_list)
all_genes_unique <- unique(all_genes)
count <- sapply(all_genes_unique, function(g) sum(sapply(sensitivity_results$top_genes_list, function(lst) g %in% lst)))

# Create an empty data frame to store gene counts
gene_counts <- data.frame(
  gene = all_genes_unique,
  count = count,
  stringsAsFactors = FALSE
)

# Sort the gene counts by count
gene_counts <- gene_counts[order(-gene_counts$count), ]


# add scores when variable = 0
counts_when_betweeness_0 <- sensitivity_results[
  sensitivity_results$betweeness_weight == 0 &
    rowSums(sensitivity_results[, c("centrality_weight", "druggability_weight", "cryptic_pocket_weight")]) != 0, ]
counts_when_betweeness_0 <- sapply(all_genes_unique, function(g) sum(sapply(counts_when_betweeness_0$top_genes_list, function(lst) g %in% lst)))
gene_counts <- cbind(gene_counts, counts_when_betweeness_0)


counts_when_centrality_0 <- sensitivity_results[
  sensitivity_results$centrality_weight == 0 &
    rowSums(sensitivity_results[, c("betweeness_weight", "druggability_weight", "cryptic_pocket_weight")]) != 0, ]
counts_when_centrality_0 <- sapply(all_genes_unique, function(g) sum(sapply(counts_when_centrality_0$top_genes_list, function(lst) g %in% lst)))
gene_counts <- cbind(gene_counts, counts_when_centrality_0)


counts_when_druggability_0 <- sensitivity_results[
  sensitivity_results$druggability_weight == 0 &
    rowSums(sensitivity_results[, c("centrality_weight", "betweeness_weight", "cryptic_pocket_weight")]) != 0, ]
counts_when_druggability_0 <- sapply(all_genes_unique, function(g) sum(sapply(counts_when_druggability_0$top_genes_list, function(lst) g %in% lst)))
gene_counts <- cbind(gene_counts, counts_when_druggability_0)


counts_when_cryptic_pockets_0 <- sensitivity_results[
  sensitivity_results$cryptic_pocket_weight == 0 &
    rowSums(sensitivity_results[, c("centrality_weight", "betweeness_weight", "druggability_weight")]) != 0, ]
counts_when_cryptic_pockets_0 <- sapply(all_genes_unique, function(g) sum(sapply(counts_when_cryptic_pockets_0$top_genes_list, function(lst) g %in% lst)))
gene_counts <- cbind(gene_counts, counts_when_cryptic_pockets_0)



# add gene description
gene_description <- getBM(attributes = c("external_gene_name", "description"), 
                          filters = "external_gene_name", 
                          values = gene_counts$gene, 
                          mart = ensembl)

final_gene_counts <- merge(gene_description, gene_counts, by.x = "external_gene_name", by.y = "gene")

final_gene_counts <- final_gene_counts[order(-final_gene_counts$count), ]

final_gene_counts$description <- gsub("\\s*\\[.*?\\]", "", final_gene_counts$description)

rownames(final_gene_counts) <- NULL


write.csv(final_gene_counts, "intermediate/paired/lumB/final_gene_counts.csv", row.names = F)
