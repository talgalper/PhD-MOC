library(TCGAbiolinks)
library(SummarizedExperiment)
library(edgeR)
library(tidyverse)
library(biomaRt)
library(PCSF)
library(progress)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")


load("RData/paired/paired_subtypes.RData")

selected_barcodes <- c(lumA_paired$normal, lumA_paired$lumA)

# drop duplciate samples
selected_barcodes <- selected_barcodes[!duplicated(selected_barcodes) & !duplicated(selected_barcodes, fromLast = TRUE)]


LumA_query <- GDCquery(project = "TCGA-BRCA",
                       access = "open", 
                       data.category = "Transcriptome Profiling",
                       experimental.strategy = "RNA-Seq",
                       barcode = selected_barcodes)

GDCdownload(LumA_query)

LumA_data <- GDCprepare(LumA_query, summarizedExperiment = T)


LumA_unstranded <- assay(LumA_data, "unstranded")  
rownames(LumA_unstranded) <- gsub("\\.\\d+", "", rownames(LumA_unstranded))
LumA_unstranded <- as.data.frame(LumA_unstranded)



subtype_subset <- master[master$cases %in% colnames(LumA_unstranded), ]
subtype_subset <- subtype_subset[match(colnames(LumA_unstranded), subtype_subset$cases), ]

group <- factor(subtype_subset$Subtype_Selected)
pat_id <- factor(subtype_subset$bcr_patient_barcode)

counts_filt <- filterByExpr(LumA_unstranded, group = group)
counts_filt <- LumA_unstranded[counts_filt, ]

low_exp_genes <- LumA_unstranded[!counts_filt, ]

#### Differential expression analysis ####

data <- DGEList(counts = counts_filt, group = group)

design <- model.matrix(~group + pat_id)

# Estimate a common negative binomial dispersion parameter for a DGE dataset with a general experimental design
common <- estimateGLMCommonDisp(data, design, verbose = T)

# Estimate the abundance-dispersion trend by Cox-Reid approximate profile likelihood.
trend <- estimateGLMTrendedDisp(common, design)

# Compute an empirical Bayes estimate of the negative binomial dispersion parameter for each tag, 
# with expression levels specified by a log-linear model.
tagwise <- estimateGLMTagwiseDisp(trend, design)

# Fit a negative binomial generalized log-linear model to the read counts for each gene. 
# Conduct genewise statistical tests for a given coefficient or coefficient contrast.
fit <- glmFit(tagwise, design)

# Conduct genewise statistical tests for a given coefficient or coefficient contrast.
lrt <- glmLRT(fit, coef = 2)

# Extract the most differentially expressed genes (or sequence tags) from a test object, 
# ranked either by p-value or by absolute log-fold-change.
toptags <- topTags(lrt, n = Inf)

# Identify which genes are significantly differentially expressed from 
# an edgeR fit object containing p-values and test statistics.
dif_exp <- decideTestsDGE(lrt, p = 0.05, adjust = "fdr", lfc = 1)
print(summary(dif_exp))

dif_exp_genes <- rownames(tagwise)[as.logical(dif_exp)]

# create a results df
hits <- toptags$table[toptags$table$FDR < 0.1, ]
hits <- rownames_to_column(hits)
rownames(hits) <- hits$rowname
colnames(hits)[1] <- "gene_id"

dif_exp <- hits[dif_exp_genes, ]

# Plot the genewise biological coefficient of variation (BCV) against gene abundance (in log2 counts per million).
plotBCV(tagwise, main = "biological coefficient of variation")

# Make a mean-difference plot of two libraries of count data with smearing of points with very low counts, 
# especially those that are zero for one of the columns.
plotSmear(lrt, de.tags = dif_exp_genes, main = "Mean-difference plot")

# plot Pvalues of different logFC scores
ggplot(hits, aes(x=logFC, y=-log(FDR))) + geom_point() + labs(title = "Adjusted logFC")

paired_hits <- toptags$table
paired_hits_subset <- paired_hits[paired_hits$FDR < 0.1, ]







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

write.table(gene_data$external_gene_name, "intermediate/paired/LumA/gene_list.txt", quote = F, row.names = F, col.names = F)


#### get interaction data ####

string_edge_data <- read.table("intermediate/paired/LumA/STRING network (physical) default edge.csv", header = T, sep = ",", stringsAsFactors = F)
ppi_list <- subset(string_edge_data, select = c("name", "stringdb..score"))
ppi_list <- ppi_list %>% 
  separate(name, sep = " ", into = c("node_1", "del", "node_2"))
ppi_list <- subset(ppi_list, select = c("node_1", "node_2", "stringdb..score"))
ppi_list$node_1 <- gsub(".*.\\.", "", ppi_list$node_1)
ppi_list$node_2 <- gsub(".*.\\.", "", ppi_list$node_2)

string_node_data <- read.table("intermediate/paired/LumA/STRING network (physical) default node.csv", header = T, sep = ",", stringsAsFactors = F)
node_list <- subset(string_node_data, select = c("name", "query.term"))
node_list$name <- gsub(".*.\\.", "", node_list$name)
ppi_list$original_order <- seq_len(nrow(ppi_list))
merged_df <- merge(ppi_list, node_list, by.x = "node_1", by.y = "name", all.x = TRUE)
merged_df <- merge(merged_df, node_list, by.x = "node_2", by.y = "name", all.x = TRUE)
merged_df <- merged_df[order(merged_df$original_order), ]

final_df <- merged_df[, c("query.term.x", "query.term.y", "stringdb..score")]
colnames(final_df) <- c("node_1", "node_2", "score")

save(final_df, gene_data, file = "RData/paired/LumA/PCSF_input.RData")


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

save(subnet, file = "RData/paired/LumA/PCSF_subnet.RData")






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

write.csv(df, "intermediate/paired/LumA/PCSF_output.csv", row.names = F)

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


# load Fpocket data
af_drugability <- read.csv("../druggability_results/fpocket_druggability.csv")

# merge PCSF data with AF
PCSF_results <- merge(PCSF_master, af_drugability, by.x = "uniprot_gn_id", by.y = "uniprot_id")

missing_genes <- unique(PCSF_master$external_gene_name)[!unique(PCSF_master$external_gene_name) %in% unique(PCSF_results$external_gene_name)]


write.csv(PCSF_results, "intermediate/paired/LumA/PCSF_druggability.csv")


# run citation scoring tool and read in
citation_scores <- read.csv("intermediate/LumA/citation_scores.csv")

PCSF_results <- merge(PCSF_results, citation_scores, by.x = "external_gene_name", by.y = "gene_id")
PCSF_results <- na.omit(PCSF_results) # for some reason the merge created empty duplicates.
PCSF_results <- PCSF_results[order(-PCSF_results$druggability), ]
rownames(PCSF_results) <- NULL


# keep structure duplicates with highest druggability
filtered_results <- PCSF_results %>%
  group_by(external_gene_name) %>%
  filter(druggability == max(druggability))

write.csv(filtered_results, "intermediate/LumA/PCSF_master_unique.csv")



#### Ranking ####
pcsf_master <- read.csv("intermediate/LumA/PCSF_master_unique.csv", row.names = 1)

pocketminer_data <- read.csv("../pocketminer/results/pocketminer_results_3.0.csv")
pcsf_master <- merge(pcsf_master, pocketminer_data, by = "ID")

missing_genes <- filtered_results$external_gene_name[!filtered_results$external_gene_name %in% pcsf_master$external_gene_name]


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


write.csv(final_gene_counts, "intermediate/paired/LumA/final_gene_counts.csv", row.names = F)
