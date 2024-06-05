library(PCSF)
library(plyr)
library(biomaRt)
library(edgeR)
library(ggplot2)
library(reshape2)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr) #need to install tidy packages separately because of biomart thing
library(tidyr)
library(tidyverse)
library(progress)
library(rDGIdb)


ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")




#### query TCGA for TCGA-OV project data ####
getProjectSummary('TCGA-BRCA')

query_TCGA <- GDCquery(project = "TCGA-BRCA",
                       access = "open", 
                       data.category = "Transcriptome Profiling",
                       experimental.strategy = "RNA-Seq")

query_output <- getResults(query_TCGA)

clinical <- GDCquery_clinic(project = "TCGA-BRCA",
                            type = "clinical")

clinical_query <- clinical[complete.cases(clinical$ajcc_pathologic_stage), ]
clinical_query <- merge(query_output, clinical_query, by.x = "cases.submitter_id", by.y = "submitter_id")
clinical_query <- subset(clinical_query, select = c("cases", "cases.submitter_id", "ajcc_pathologic_stage", 
                                                    "tissue_or_organ_of_origin", "sample_type"))

table(clinical_query$ajcc_pathologic_stage)



subtypes <- PanCancerAtlas_subtypes()

#common <- query_output[query_output$cases %in% subtypes$pan.samplesID, ]

common <- merge(clinical_query, subtypes, by.x = "cases", by.y = "pan.samplesID")
common <- subset(common, select = c("cases", "Subtype_Selected", "sample_type", "ajcc_pathologic_stage"))


selected_barcodes <- common[common$Subtype_Selected == "BRCA.LumA", ]

LumA_query <- GDCquery(project = "TCGA-BRCA",
                        access = "open", 
                        data.category = "Transcriptome Profiling",
                        experimental.strategy = "RNA-Seq",
                        barcode = selected_barcodes$cases)

GDCdownload(LumA_query)

LumA_data <- GDCprepare(LumA_query, summarizedExperiment = T)


LumA_unstranded <- assay(LumA_data, "unstranded")  
rownames(LumA_unstranded) <- gsub("\\.\\d+", "", rownames(LumA_unstranded))
LumA_unstranded <- as.data.frame(LumA_unstranded)



barplot(colSums(LumA_unstranded),
        xaxt = "n")

barplot(colSums(normal_unstranded),
        xaxt = "n")


#### Remove low activity genes ####
colnames(LumA_unstranded) <- paste("LumA", 1:ncol(LumA_unstranded), sep = "_")
colnames(normal_unstranded) <- paste("normal", 1:ncol(normal_unstranded), sep = "_")

merged_df <- merge(LumA_unstranded, normal_unstranded, by = "row.names")

merged_df <- column_to_rownames(merged_df, "Row.names")

gene_id <- rownames(merged_df)

# subset genes that have a cpm of greater than one, in more than 50% of the samples
keepTheseGenes <- (rowSums(cpm(merged_df) > 1) >= ncol(merged_df)/2)
print(summary(keepTheseGenes))

# add gene ids back into df
merged_df <- cbind(gene_id, merged_df)

removedGenes <- merged_df$gene_id[!keepTheseGenes]
removedGenes <- as.data.frame(removedGenes)
colnames(removedGenes)[1] <- "gene_id"

merged_df <- merged_df[keepTheseGenes, ]

merged_df <- merged_df[, -1]


# top x% of samples
#num_genes <- nrow(merged_df)
#top_threshold <- ceiling(num_genes * 0.50)
#
## min required samples
#num_samples <- ncol(merged_df)
#min_required_samples <- ceiling(num_samples * 0.10)
#
#
#
## Extract the gene names and data columns
#gene_names <- rownames(merged_df)
#data_columns <- merged_df[, 2:ncol(merged_df)]
#
## Initialise an empty data frame
#results_df <- data.frame(gene_id = gene_names)
#
## Iterate through each data column and perform the test
#for (col_id in seq_along(data_columns)) {
#  # Get the column name and the data
#  col_name <- colnames(data_columns)[col_id]
#  col_data <- data_columns[, col_id]
#  
#  # Sort the data and mark genes within the top_threshold rows as TRUE
#  sorted_indices <- order(col_data, decreasing = TRUE)
#  top_indices <- sorted_indices[1:top_threshold]
#  
#  # Create a logical vector for gene presence in top_threshold
#  gene_presence <- rep(FALSE, nrow(merged_df))
#  gene_presence[top_indices] <- TRUE
#  
#  # Add to the results data frame
#  results_df <- cbind(results_df, gene_presence)
#}
#
## Rename columns 
#colnames(results_df) <- c("gene_id", colnames(data_columns))
#
#
## Calculate the number of TRUE values for each gene across columns
#gene_counts <- rowSums(results_df[, -1]) # Exclude the first column ("gene_id")
#failed_genes <- gene_counts < min_required_samples
#print(summary(failed_genes))
#
#BRCA_data <- subset(merged_df, !(rownames(merged_df) %in% gene_names[failed_genes]))







#### Differential expression analysis ####

# Select the columns that start with "LumA"
cols <- grep("^LumA", colnames(merged_df))

# Assign the unstranded columns to the cancer group
cancer <- colnames(merged_df)[cols]

# Assign the rest of the columns to the healty group
benign <- setdiff(colnames(merged_df), cancer)

# Create the group variable
group <- factor(c(rep("cancer", length(cancer)), rep("benign", length(benign))))

data <- DGEList(counts = merged_df, group = group)

design <- model.matrix(~group)

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
dif_exp <- decideTests(lrt, p = 0.05, adjust = "fdr", lfc = 1)
print(summary(dif_exp))

dif_exp_genes <- rownames(tagwise)[as.logical(dif_exp)]

# create a results df
hits <- toptags$table[toptags$table$FDR < 0.1, ]
colnames <- colnames(hits)
hits$gene_id <- rownames(hits)
hits <- hits[,c("gene_id", colnames)]

dif_exp <- hits[dif_exp_genes, ]

write.csv(hits, "intermediate/LumA/DE_results.csv", row.names = F)

# Plot the genewise biological coefficient of variation (BCV) against gene abundance (in log2 counts per million).
plotBCV(tagwise, main = "biological coefficient of variation")

# Make a mean-difference plot of two libraries of count data with smearing of points with very low counts, 
# especially those that are zero for one of the columns.
plotSmear(lrt, de.tags = dif_exp_genes, main = "Mean-difference plot")

# plot Pvalues of different logFC scores
ggplot(hits, aes(x=logFC, y=-log(FDR))) + geom_point() + labs(title = "Adjusted logFC")



#### Convert to Gene IDs and get interaction data ####

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

write.table(gene_data$external_gene_name, "intermediate/LumA/gene_list.txt", quote = F, row.names = F, col.names = F)


#### get interaction data ####

string_edge_data <- read.table("intermediate/LumA/STRING network (physical) default edge.csv", header = T, sep = ",", stringsAsFactors = F)
ppi_list <- subset(string_edge_data, select = c("name", "stringdb..score"))
ppi_list <- ppi_list %>% 
  separate(name, sep = " ", into = c("node_1", "del", "node_2"))
ppi_list <- subset(ppi_list, select = c("node_1", "node_2", "stringdb..score"))
ppi_list$node_1 <- gsub(".*.\\.", "", ppi_list$node_1)
ppi_list$node_2 <- gsub(".*.\\.", "", ppi_list$node_2)

string_node_data <- read.table("intermediate/LumA/STRING network (physical) default node.csv", header = T, sep = ",", stringsAsFactors = F)
node_list <- subset(string_node_data, select = c("name", "query.term"))
node_list$name <- gsub(".*.\\.", "", node_list$name)
ppi_list$original_order <- seq_len(nrow(ppi_list))
merged_df <- merge(ppi_list, node_list, by.x = "node_1", by.y = "name", all.x = TRUE)
merged_df <- merge(merged_df, node_list, by.x = "node_2", by.y = "name", all.x = TRUE)
merged_df <- merged_df[order(merged_df$original_order), ]

final_df <- merged_df[, c("query.term.x", "query.term.y", "stringdb..score")]
colnames(final_df) <- c("node_1", "node_2", "score")

save(final_df, gene_data, file = "RData/LumA/PCSF_input.RData")


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

save(subnet, file = "RData/LumA/PCSF_subnet.RData")






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

write.csv(df, "intermediate/LumA/PCSF_output.csv", row.names = F)

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


write.csv(PCSF_results, "intermediate/LumA/PCSF_druggability.csv")


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


write.csv(final_gene_counts, "intermediate/LumA/final_gene_counts.csv", row.names = F)





