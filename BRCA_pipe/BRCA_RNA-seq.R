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
library(GEOquery)
library(RobustRankAggreg)
library(progress)
library(rDGIdb)


ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")




#### query TCGA for TCGA-OV project data ####
getProjectSummary('TCGA-BRCA')

query_TCGA <- GDCquery(project = "TCGA-BRCA",
                       access = "open", 
                       data.category = "Transcriptome Profiling")

query_output <- getResults(query_TCGA)

clinical <- GDCquery_clinic(project = "TCGA-BRCA",
                            type = "clinical")

clinical <- merge(query_output, clinical, by.x = "cases.submitter_id", by.y = "submitter_id")
clinical <- subset(clinical, select = c("cases.submitter_id", "ajcc_pathologic_stage", "tissue_or_organ_of_origin"))

table(clinical$ajcc_pathologic_stage)



subtypes <- PanCancerAtlas_subtypes()

common <- query_output[query_output$cases %in% subtypes$pan.samplesID, ]

common <- merge(query_output, subtypes, by.x = "cases", by.y = "pan.samplesID")
selected_barcodes <- subset(common, select = c("cases", "Subtype_Selected", "sample_type"))
selected_barcodes <- selected_barcodes[selected_barcodes$Subtype_Selected == "BRCA.LumB", ]


LumB_query <- GDCquery(project = "TCGA-BRCA",
                       access = "open", 
                       data.category = "Transcriptome Profiling",
                       barcode = selected_barcodes$cases)

GDCdownload(LumB_query)

LumB_data <- GDCprepare(LumB_query, summarizedExperiment = T)


LumB_unstranded <- assay(LumB_data, "unstranded")  
rownames(LumB_unstranded) <- gsub("\\.\\d+", "", rownames(LumB_unstranded))
LumB_unstranded <- as.data.frame(LumB_unstranded)



GTEx_raw <- read.table("bulk-gex_v8_rna-seq_counts-by-tissue_gene_reads_2017-06-05_v8_breast_mammary_tissue.gct", skip = 2, header = T)

GTEx_data <- GTEx_raw[, -c(1:3)]
rownames(GTEx_data) <- GTEx_raw$Name
rownames(GTEx_data) <- gsub("\\.\\d+", "", rownames(GTEx_data))


save(LumB_unstranded, GTEx_data, file = "RData/BRCA_DE_data.RData")

barplot(colSums(LumB_unstranded),
        xaxt = "n")

barplot(colSums(GTEx_data),
        xaxt = "n")

#### Remove low activity genes ####
merged_df <- merge(LumB_unstranded, GTEx_data, by = "row.names")

merged_df <- column_to_rownames(merged_df, "Row.names")

gene_id <- rownames(merged_df)

# subset genes that have a cpm of greater than one, in more than 50% of the samples
#keepTheseGenes <- (rowSums(cpm(merged_df) > 1) >= ncol(merged_df)/2)
#print(summary(keepTheseGenes))



# top x% of samples
num_genes <- nrow(merged_df)
top_threshold <- ceiling(num_genes * 0.50)

# min required samples
num_samples <- ncol(merged_df)
min_required_samples <- ceiling(num_samples * 0.10)



# Extract the gene names and data columns
gene_names <- rownames(merged_df)
data_columns <- merged_df[, 2:ncol(merged_df)]

# Initialise an empty data frame
results_df <- data.frame(gene_id = gene_names)

# Iterate through each data column and perform the test
for (col_id in seq_along(data_columns)) {
  # Get the column name and the data
  col_name <- colnames(data_columns)[col_id]
  col_data <- data_columns[, col_id]
  
  # Sort the data and mark genes within the top_threshold rows as TRUE
  sorted_indices <- order(col_data, decreasing = TRUE)
  top_indices <- sorted_indices[1:top_threshold]
  
  # Create a logical vector for gene presence in top_threshold
  gene_presence <- rep(FALSE, nrow(merged_df))
  gene_presence[top_indices] <- TRUE
  
  # Add to the results data frame
  results_df <- cbind(results_df, gene_presence)
}

# Rename columns 
colnames(results_df) <- c("gene_id", colnames(data_columns))


# Calculate the number of TRUE values for each gene across columns
gene_counts <- rowSums(results_df[, -1]) # Exclude the first column ("gene_id")
failed_genes <- gene_counts < min_required_samples
print(summary(failed_genes))

BRCA_data <- subset(merged_df, !(rownames(merged_df) %in% gene_names[failed_genes]))


# add gene ids back into df
#merged_df <- cbind(gene_id, merged_df)
#
#removedGenes <- merged_df$gene_id[!keepTheseGenes]
#removedGenes <- as.data.frame(removedGenes)
#colnames(removedGenes)[1] <- "gene_id"
#
#merged_df <- merged_df[keepTheseGenes, ]
#
#merged_df <- merged_df[, -1]




#### Differential expression analysis ####

# Select the columns that start with "TCGA"
cols <- grep("^TCGA", colnames(BRCA_data))

# Assign the unstranded columns to the cancer group
cancer <- colnames(BRCA_data)[cols]

# Assign the rest of the columns to the healty group
benign <- setdiff(colnames(BRCA_data), cancer)

# Create the group variable
group <- factor(c(rep("cancer", length(cancer)), rep("benign", length(benign))))

data <- DGEList(counts = BRCA_data, group = group)

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
dif_exp <- decideTestsDGE(lrt, p = 0.05, adjust = "fdr", lfc = 2)
print(summary(dif_exp))

dif_exp_genes <- rownames(tagwise)[as.logical(dif_exp)]

# create a results df
hits <- toptags$table[toptags$table$FDR < 0.1, ]
colnames <- colnames(hits)
hits$gene_id <- rownames(hits)
hits <- hits[,c("gene_id", colnames)]

dif_exp <- hits[dif_exp_genes, ]


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

write.table(gene_data$external_gene_name, "intermediate/gene_list.txt", quote = F, row.names = F, col.names = F)


#### get interaction data ####

string_edge_data <- read.table("intermediate/STRING network (physical) default edge.csv", header = T, sep = ",", stringsAsFactors = F)
ppi_list <- subset(string_edge_data, select = c("name", "stringdb..score"))
ppi_list <- ppi_list %>% 
  separate(name, sep = " ", into = c("node_1", "del", "node_2"))
ppi_list <- subset(ppi_list, select = c("node_1", "node_2", "stringdb..score"))
ppi_list$node_1 <- gsub(".*.\\.", "", ppi_list$node_1)
ppi_list$node_2 <- gsub(".*.\\.", "", ppi_list$node_2)

string_node_data <- read.table("intermediate/STRING network (physical) default node.csv", header = T, sep = ",", stringsAsFactors = F)
node_list <- subset(string_node_data, select = c("name", "query.term"))
node_list$name <- gsub(".*.\\.", "", node_list$name)
ppi_list$original_order <- seq_len(nrow(ppi_list))
merged_df <- merge(ppi_list, node_list, by.x = "node_1", by.y = "name", all.x = TRUE)
merged_df <- merge(merged_df, node_list, by.x = "node_2", by.y = "name", all.x = TRUE)
merged_df <- merged_df[order(merged_df$original_order), ]

final_df <- merged_df[, c("query.term.x", "query.term.y", "stringdb..score")]
colnames(final_df) <- c("node_1", "node_2", "score")

save(final_df, gene_data, file = "RData/PCSF_input.RData")


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

save(subnet, file = "RData/PCSF_subnet.RData")






# extract cluster data
clust <- clusters(subnet)
df <- data.frame(gene_id = names(clust$membership), cluster = factor(clust$membership))
betweenness <- betweenness(subnet) 
centrality <- degree(subnet) 
df$betweenness <- betweenness[as.character(df$gene_id)]
df$degree_centrality <- centrality[as.character(df$gene_id)]
df$betweenness <- as.integer(df$betweenness)
df$degree_centrality <- as.integer(df$degree_centrality)

rownames(df) <- 1:nrow(df)

df <- df[order(-df$degree_centrality), ]


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

write.csv(PCSF_results, "intermediate/PCSF_druggability.csv")


# run citation scoring tool and read in
citation_scores <- read.csv("intermediate/citation_scores.csv")

PCSF_results <- merge(PCSF_results, citation_scores, by.x = "external_gene_name", by.y = "gene_id")
PCSF_results <- na.omit(PCSF_results) # for some reason the merge created empty duplicates.
PCSF_results <- PCSF_results[order(-PCSF_results$druggability), ]
rownames(PCSF_results) <- NULL


# Jess wants to keep structure duplicates with highest druggability
filtered_results <- PCSF_results %>%
  group_by(external_gene_name) %>%
  filter(druggability == max(druggability))

write.csv(filtered_results, "intermediate/PCSF_master_unique.csv")



#### PocketMiner data ####
pocketminer_data <- read.csv("../pocketminer/results/pocketminer_results_3.0.csv")

















#### Druggability rank ####

data <- read.csv("intermediate/PCSF_master_unique.csv", row.names = 1)
rownames(data) <- NULL
master <- merge(data, pocketminer_data, by = "ID")
master <- subset(master, select = c("external_gene_name", "druggability", "struct_score", "max_hit"))

# normalise pocketminer scores
pocketminer_norm <- (master$max_hit - min(master$max_hit)) / 
  (max(master$max_hit) - min(master$max_hit))

# convert structure scores to decimal values from percentages
structure_norm <- master$struct_score / 100



# Define the number of weight values (genes in top x) to test (including 0)
# i.e. 0, 0.1, 0.2, 0.3, ..., 1
num_weights <- 11

# Create weight value sequences that sum up to 1
fpocket_w_values <- seq(0, 1, length.out = num_weights)
pocketminer_w_values <- seq(0, 1, length.out = num_weights)
structure_w_values <- seq(0, 1, length.out = num_weights)


# Calculate the total number of iterations
total_iterations <- num_weights ^ 3

# Create a progress bar
pb <- progress_bar$new(total = total_iterations)

# Create an empty data frame with proper column names
sensitivity_results <- data.frame(
  fpocket_weight = numeric(),
  pocketminer_weight = numeric(),
  structure_weight = numeric(),
  top_genes = character(),
  stringsAsFactors = FALSE
)

# Perform sensitivity test
for (fpocket_w in fpocket_w_values) {
  for (pocketminer_w in pocketminer_w_values) {
    for (structure_w in structure_w_values) {
      # Check if the weights sum up to 1
      if (fpocket_w + pocketminer_w + structure_w == 1) {
        # Combine scores using current weights
        combined_score <- fpocket_w * master$druggability +
          pocketminer_w * pocketminer_norm +
          structure_w * structure_norm
        
        # Rank the genes based on the combined score
        master_ranked <- master
        master_ranked$combined_score <- combined_score
        master_ranked <- master_ranked[order(-master_ranked$combined_score), ]
        
        # Select the top 10 genes
        top_genes <- head(master_ranked$external_gene_name, 10)
        
        # Add results to the sensitivity_results data frame
        sensitivity_results <- rbind(sensitivity_results, list(
          fpocket_weight = fpocket_w,
          pocketminer_weight = pocketminer_w,
          structure_weight = structure_w,
          top_genes = paste(top_genes, collapse = ', ')
        ))
      }
      
      # Increment the progress bar
      pb$tick()
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

counts_when_fpocket_0 <- sensitivity_results[
  sensitivity_results$fpocket_weight == 0 &
    rowSums(sensitivity_results[, c("fpocket_weight", "pocketminer_weight", "structure_weight")]) != 0, 
]
counts_when_fpocket_0 <- sapply(all_genes_unique, function(g) sum(sapply(counts_when_fpocket_0$top_genes_list, function(lst) g %in% lst)))

gene_counts <- cbind(gene_counts, counts_when_fpocket_0)

counts_when_pocketminer_0 <- sensitivity_results[
  sensitivity_results$pocketminer_weight == 0 &
    rowSums(sensitivity_results[, c("fpocket_weight", "fpocket_weight", "structure_weight")]) != 0, 
]
counts_when_pocketminer_0 <- sapply(all_genes_unique, function(g) sum(sapply(counts_when_pocketminer_0$top_genes_list, function(lst) g %in% lst)))

gene_counts <- cbind(gene_counts, counts_when_pocketminer_0)

counts_when_structure_0 <- sensitivity_results[
  sensitivity_results$structure_weight == 0 &
    rowSums(sensitivity_results[, c("fpocket_weight", "pocketminer_weight", "structure_weight")]) != 0, 
]
counts_when_structure_0 <- sapply(all_genes_unique, function(g) sum(sapply(counts_when_structure_0$top_genes_list, function(lst) g %in% lst)))

final_gene_counts <- cbind(gene_counts, counts_when_structure_0)



# add gene description
gene_description <- getBM(attributes = c("external_gene_name", "description"), 
                          filters = "external_gene_name", 
                          values = final_gene_counts$gene, 
                          mart = ensembl)

final_gene_counts <- merge(gene_description, final_gene_counts, by.x = "external_gene_name", by.y = "gene")

final_gene_counts$description <- gsub("\\s*\\[.*?\\]", "", final_gene_counts$description)


DGIdb <- queryDGIdb(final_gene_counts$external_gene_name)
results <- byGene(DGIdb)
results <- subset(results, select = c("Gene", "DistinctDrugCount"))
final_gene_counts <- merge(final_gene_counts, results, by.x = "external_gene_name", by.y = "Gene")
final_gene_counts <- final_gene_counts[order(-final_gene_counts$count), ]
rownames(final_gene_counts) <- NULL








cit_scores <- read.csv("intermediate/citation_scores.csv")
master <- merge(master, cit_scores, by.x = "external_gene_name", by.y = "gene_id")

citation_norm <- log(master$MeSH_count + 1)
citation_norm <- (citation_norm - min(citation_norm)) / (max(citation_norm) - min(citation_norm))

# Create weight value sequences that sum up to 1
fpocket_w_values <- seq(0, 1, length.out = num_weights)
pocketminer_w_values <- seq(0, 1, length.out = num_weights)
structure_w_values <- seq(0, 1, length.out = num_weights)
citation_w_values <- seq(0, 1, length.out = num_weights)


# Calculate the total number of iterations
total_iterations <- num_weights ^ 4

# Create a progress bar
pb <- progress_bar$new(total = total_iterations)

# Create an empty data frame with proper column names
sensitivity_results <- data.frame(
  fpocket_weight = numeric(),
  pocketminer_weight = numeric(),
  structure_weight = numeric(),
  citation_weight = numeric(),
  top_genes = character(),
  stringsAsFactors = FALSE
)

# Perform sensitivity test
for (fpocket_w in fpocket_w_values) {
  for (pocketminer_w in pocketminer_w_values) {
    for (structure_w in structure_w_values) {
      for (citation_w in citation_w_values) {
        # Check if the weights sum up to 1
        if (fpocket_w + pocketminer_w + structure_w + citation_w == 1) {
          # Combine scores using current weights
          combined_score <- fpocket_w * master$druggability +
            pocketminer_w * pocketminer_norm +
            structure_w * structure_norm -
            citation_w * citation_norm
          
          
          # Rank the genes based on the combined score
          master_ranked <- master
          master_ranked$combined_score <- combined_score
          master_ranked <- master_ranked[order(-master_ranked$combined_score), ]
          
          # Select the top 10 genes
          top_genes <- head(master_ranked$external_gene_name, 10)
          
          # Add results to the sensitivity_results data frame
          sensitivity_results <- rbind(sensitivity_results, list(
            fpocket_weight = fpocket_w,
            pocketminer_weight = pocketminer_w,
            structure_weight = structure_w,
            citation_weight = citation_w,
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



counts_when_fpocket_0 <- sensitivity_results[
  sensitivity_results$fpocket_weight == 0 &
    rowSums(sensitivity_results[, c("citation_weight", "pocketminer_weight", "structure_weight")]) != 0, 
]
counts_when_fpocket_0 <- sapply(all_genes_unique, function(g) sum(sapply(counts_when_fpocket_0$top_genes_list, function(lst) g %in% lst)))

gene_counts <- cbind(gene_counts, counts_when_fpocket_0)


counts_when_pocketminer_0 <- sensitivity_results[
  sensitivity_results$pocketminer_weight == 0 &
    rowSums(sensitivity_results[, c("fpocket_weight", "citation_weight", "structure_weight")]) != 0, 
]
counts_when_pocketminer_0 <- sapply(all_genes_unique, function(g) sum(sapply(counts_when_pocketminer_0$top_genes_list, function(lst) g %in% lst)))

gene_counts <- cbind(gene_counts, counts_when_pocketminer_0)


counts_when_structure_0 <- sensitivity_results[
  sensitivity_results$structure_weight == 0 &
    rowSums(sensitivity_results[, c("fpocket_weight", "pocketminer_weight", "citation_weight")]) != 0, 
]
counts_when_structure_0 <- sapply(all_genes_unique, function(g) sum(sapply(counts_when_structure_0$top_genes_list, function(lst) g %in% lst)))

gene_counts <- cbind(gene_counts, counts_when_structure_0)


counts_when_citation_0 <- sensitivity_results[
  sensitivity_results$citation_weight == 0 &
    rowSums(sensitivity_results[, c("fpocket_weight", "pocketminer_weight", "structure_weight")]) != 0, 
]
counts_when_citation_0 <- sapply(all_genes_unique, function(g) sum(sapply(counts_when_citation_0$top_genes_list, function(lst) g %in% lst)))

gene_counts_citation <- cbind(gene_counts, counts_when_citation_0)


# add gene description
gene_description <- getBM(attributes = c("external_gene_name", "description"), 
                          filters = "external_gene_name", 
                          values = gene_counts_citation$gene, 
                          mart = ensembl)

gene_counts_citation <- merge(gene_description, gene_counts_citation, by.x = "external_gene_name", by.y = "gene")

gene_counts_citation$description <- gsub("\\s*\\[.*?\\]", "", gene_counts_citation$description)


DGIdb <- queryDGIdb(gene_counts_citation$external_gene_name)
results <- byGene(DGIdb)
results <- subset(results, select = c("Gene", "DistinctDrugCount"))
gene_counts_citation <- merge(gene_counts_citation, results, by.x = "external_gene_name", by.y = "Gene")
gene_counts_citation <- gene_counts_citation[order(-gene_counts_citation$count), ]
rownames(gene_counts_citation) <- NULL


save(final_gene_counts, gene_counts_citation, file = "RData/druggability_rank.RData")




#### normalise counts for druggability score ####

norm_scores <- data.frame(gene = final_gene_counts$external_gene_name,
                          description = final_gene_counts$description,
                          counts_norm = (final_gene_counts$count - min(final_gene_counts$count)) / (max(final_gene_counts$count) - min(final_gene_counts$count)),
                          DistinctDrugCount = final_gene_counts$DistinctDrugCount)


norm_scores_cit <- data.frame(gene = gene_counts_citation$external_gene_name,
                              description = gene_counts_citation$description,
                              counts_norm = (gene_counts_citation$count - min(gene_counts_citation$count)) / (max(gene_counts_citation$count) - min(gene_counts_citation$count)),
                              DistinctDrugCount = gene_counts_citation$DistinctDrugCount)




#### Rank sensitivity ####
PCSF_master <- read.csv("intermediate/PCSF_master_unique.csv", row.names = 1)
PCSF_master <- merge(PCSF_master, pocketminer_data, by = "ID")

betweeness_norm <- (PCSF_master$betweenness - min(PCSF_master$betweenness)) / (max(PCSF_master$betweenness) - min(PCSF_master$betweenness))

centrality_norm <- (PCSF_master$degree_centrality - min(PCSF_master$degree_centrality)) / (max(PCSF_master$degree_centrality) - min(PCSF_master$degree_centrality))

# log transformation
citation_norm <- log(PCSF_master$MeSH_count + 1)
citation_norm <- (citation_norm - min(citation_norm)) / (max(citation_norm) - min(citation_norm))




# Define the number of weight values (genes in top x) to test (including 0)
# i.e. 0, 0.1, 0.2, 0.3, ..., 1
num_weights <- 11

# Create weight value sequences that sum up to 1
betweeness_w_values <- seq(0, 1, length.out = num_weights)
citation_w_values <- seq(0, 1, length.out = num_weights)
centrality_w_values <- seq(0, 1, length.out = num_weights)
druggability_w_values <- seq(0, 1, length.out = num_weights)
cryptic_pocket_w_values <- seq(0, 1, length.out = num_weights)


# Calculate the total number of iterations a.k.a number of variables
total_iterations <- num_weights ^ 5

# Create a progress bar
pb <- progress_bar$new(total = total_iterations)

# Create an empty data frame with proper column names
sensitivity_results <- data.frame(
  betweeness_weight = numeric(),
  citation_weight = numeric(),
  centrality_weight = numeric(),
  druggability_weight = numeric(),
  cryptic_pocket_weight = numeric(),
  top_genes = character(),
  stringsAsFactors = FALSE
)

# Perform sensitivity test
for (betweeness_w in betweeness_w_values) {
  for (citation_w in citation_w_values) {
    for (centrality_w in centrality_w_values) {
      for (druggability_w in druggability_w_values) {
        for (cryptic_pocket_w in cryptic_pocket_w_values) {
          # Check if the weights sum up to 1
          if (betweeness_w + citation_w + centrality_w + druggability_w + cryptic_pocket_w  == 1) {
            # Combine scores using current weights
            combined_score <- (betweeness_w * betweeness_norm) +
              (centrality_w * centrality_norm) +
              (druggability_w * PCSF_master$druggability) +
              (cryptic_pocket_w * PCSF_master$max_hit) -
              (citation_w * citation_norm)
            
            # Rank the genes based on the combined score
            PCSF_master_ranked <- PCSF_master
            PCSF_master_ranked$combined_score <- combined_score
            PCSF_master_ranked <- PCSF_master_ranked[order(-PCSF_master_ranked$combined_score), ]
            
            # Select the top 10 genes
            top_genes <- head(PCSF_master_ranked$external_gene_name, 10)
            
            # Add results to the sensitivity_results data frame
            sensitivity_results <- rbind(sensitivity_results, list(
              betweeness_weight = betweeness_w,
              citation_weight = citation_w,
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
    rowSums(sensitivity_results[, c("citation_weight", "centrality_weight", "druggability_weight", "cryptic_pocket_weight")]) != 0, ]
counts_when_betweeness_0 <- sapply(all_genes_unique, function(g) sum(sapply(counts_when_betweeness_0$top_genes_list, function(lst) g %in% lst)))
gene_counts <- cbind(gene_counts, counts_when_betweeness_0)


counts_when_centrality_0 <- sensitivity_results[
  sensitivity_results$centrality_weight == 0 &
    rowSums(sensitivity_results[, c("citation_weight", "betweeness_weight", "druggability_weight", "cryptic_pocket_weight")]) != 0, ]
counts_when_centrality_0 <- sapply(all_genes_unique, function(g) sum(sapply(counts_when_centrality_0$top_genes_list, function(lst) g %in% lst)))
gene_counts <- cbind(gene_counts, counts_when_centrality_0)


counts_when_citation_0 <- sensitivity_results[
  sensitivity_results$citation_weight == 0 &
    rowSums(sensitivity_results[, c("betweeness_weight", "centrality_weight", "druggability_weight", "cryptic_pocket_weight")]) != 0, ]
counts_when_citation_0 <- sapply(all_genes_unique, function(g) sum(sapply(counts_when_citation_0$top_genes_list, function(lst) g %in% lst)))
gene_counts <- cbind(gene_counts, counts_when_citation_0)


counts_when_druggability_0 <- sensitivity_results[
  sensitivity_results$druggability_weight == 0 &
    rowSums(sensitivity_results[, c("citation_weight", "centrality_weight", "betweeness_weight", "cryptic_pocket_weight")]) != 0, ]
counts_when_druggability_0 <- sapply(all_genes_unique, function(g) sum(sapply(counts_when_druggability_0$top_genes_list, function(lst) g %in% lst)))
gene_counts <- cbind(gene_counts, counts_when_druggability_0)


counts_when_cryptic_pockets_0 <- sensitivity_results[
  sensitivity_results$cryptic_pocket_weight == 0 &
    rowSums(sensitivity_results[, c("citation_weight", "centrality_weight", "betweeness_weight", "druggability_weight")]) != 0, ]
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

DGIdb <- queryDGIdb(final_gene_counts$external_gene_name)
results <- byGene(DGIdb)
results <- subset(results, select = c("Gene", "DistinctDrugCount"))
final_gene_counts <- merge(final_gene_counts, results, by.x = "external_gene_name", by.y = "Gene")
final_gene_counts <- final_gene_counts[order(-final_gene_counts$count), ]
rownames(final_gene_counts) <- NULL
























# Robust Rank Agrregation
betweenness <- subset(filtered_results, select = c("ID", "betweenness"))
betweenness <- distinct(betweenness)
betweenness$betweenness <- (betweenness$betweenness - min(betweenness$betweenness)) / (max(betweenness$betweenness) - min(betweenness$betweenness))
betweenness <- betweenness[order(betweenness$betweenness, decreasing = T), ]
betweenness <- betweenness$ID

centrality <- subset(filtered_results, select = c("ID", "degree_centrality"))
centrality <- distinct(centrality)
centrality$degree_centrality <- (centrality$degree_centrality - min(centrality$degree_centrality)) / (max(centrality$degree_centrality) - min(centrality$degree_centrality))
centrality <- centrality[order(centrality$degree_centrality, decreasing = T), ]
centrality <- centrality$ID

# log transformation of citation scores
citation <- subset(filtered_results, select = c("ID", "citation_score"))
citation <- distinct(citation)
citation$citation_score <- log(citation$citation_score + 1)
citation$citation_score <- (citation$citation_score - min(citation$citation_score)) / (max(citation$citation_score) - min(citation$citation_score))
citation <- citation[order(citation$citation_score, decreasing = F), ]
citation <- citation$ID

druggability_data <- subset(filtered_results, select = c("ID", "druggability", "num_drug_pockets"))
druggability_data <- distinct(druggability_data)
druggability <- druggability_data[order(druggability_data$druggability, decreasing = T), ]
druggability <- druggability$ID

num_drug_pockets <- druggability_data[order(druggability_data$num_drug_pockets, decreasing = T), ]
num_drug_pockets <- num_drug_pockets$ID

# add the pocketminer data as another ranking
pocketminer_data <- read.csv("../pocketminer/results/pocketminer_results.csv")
pocketminer_data <- pocketminer_data[pocketminer_data$ID %in% filtered_results$ID, ]

#pocketminer_gn_id <- getBM(attributes = c("uniprot_gn_id", "external_gene_name"), 
#                           filters = "uniprot_gn_id", 
#                           values = pocketminer_data$uniprot_id, 
#                           mart = ensembl)

#pocketminer_data <- merge(pocketminer_gn_id, pocketminer_data, by.x = "uniprot_gn_id", by.y = "uniprot_id")
#pocketminer_data[pocketminer_data == ""] <- NA
#pocketminer_data <- na.omit(pocketminer_data)
pocketminer_data <- pocketminer_data[order(pocketminer_data$max_hit, decreasing = T), ]
cryptic_pockets <- pocketminer_data$ID

# reorder for num_hits
pocketminer_data <- pocketminer_data[order(pocketminer_data$num_hits, decreasing = T), ]
num_cryp_pockets <- pocketminer_data$ID

rankings <- list(betweenness, centrality, druggability, num_drug_pockets, cryptic_pockets, num_cryp_pockets)


aggregate_ranks <- aggregateRanks(glist = rankings)

# merge gene names back in
temp <- subset(filtered_results, select = c("external_gene_name", "ID"))
aggregate_ranks <- merge(temp, aggregate_ranks, by.x = "ID", by.y = "Name")
rm(temp)

aggregate_ranks <- aggregate_ranks[order(aggregate_ranks$Score), ]
rownames(aggregate_ranks) <- NULL

# remove scientific notation
options(scipen=999)

