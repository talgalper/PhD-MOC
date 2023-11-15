library(PCSF)
library(plyr)
library(biomaRt)
library(edgeR)
library(ggplot2)
library(reshape2)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(tidyverse)
library(GEOquery)
library(RobustRankAggreg)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")




# query TCGA for TCGA-OV project data
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



merged_df <- merge(LumB_unstranded, GTEx_data, by = "row.names")

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






# Select the columns that start with "TCGA"
cols <- grep("^TCGA", colnames(merged_df))

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






data <- subset(dif_exp, select = c("gene_id", "logFC"))

# convert to gene symbol
uniprot_id <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), 
                    filters = "ensembl_gene_id", 
                    values = data$gene_id, 
                    mart = ensembl)

# remove empty rows
uniprot_id <- subset(uniprot_id, external_gene_name != "")

# check for duplicate uniprot ids
uniprot_id <- distinct(uniprot_id)

# merge back with original data
colnames(uniprot_id)[1] <- "gene_id"

uniprot_data <- merge(uniprot_id, data, by = "gene_id")

uniprot_data <- subset(uniprot_data, select = c("external_gene_name", "logFC"))

# check to see if genes that at all genes were converted at least once
missing_genes <- anti_join(data, uniprot_id, by = "gene_id")

#write.table(uniprot_data$external_gene_name, "~/Desktop/gene_list.txt", quote = F, row.names = F, col.names = F)




string_edge_data <- read.table("STRING network (physical) default edge.csv", header = T, sep = ",", stringsAsFactors = F)
ppi_list <- subset(string_edge_data, select = c("name", "stringdb..score"))
ppi_list <- ppi_list %>% 
  separate(name, sep = " ", into = c("node_1", "del", "node_2"))
ppi_list <- subset(ppi_list, select = c("node_1", "node_2", "stringdb..score"))
ppi_list$node_1 <- gsub(".*.\\.", "", ppi_list$node_1)
ppi_list$node_2 <- gsub(".*.\\.", "", ppi_list$node_2)

string_node_data <- read.table("STRING network (physical) default node.csv", header = T, sep = ",", stringsAsFactors = F)
node_list <- subset(string_node_data, select = c("name", "query.term"))
node_list$name <- gsub(".*.\\.", "", node_list$name)
ppi_list$original_order <- seq_len(nrow(ppi_list))
merged_df <- merge(ppi_list, node_list, by.x = "node_1", by.y = "name", all.x = TRUE)
merged_df <- merge(merged_df, node_list, by.x = "node_2", by.y = "name", all.x = TRUE)
merged_df <- merged_df[order(merged_df$original_order), ]

final_df <- merged_df[, c("query.term.x", "query.term.y", "stringdb..score")]
colnames(final_df) <- c("node_1", "node_2", "score")

#write.csv(final_df, "no_gepia/interaction_score.csv")


# create gene_id uniprot_id ref list.
ref_list <- getBM(attributes = c("external_gene_name","uniprot_gn_id"),
                  filters = "external_gene_name", 
                  values = node_list$query.term, 
                  mart = ensembl)


save(final_df, uniprot_data, file = "RData/PCSF_input.RData")


# set seed for reproducibility 
set.seed(1234)

# construct interactome
ppi <- construct_interactome(final_df)

# set terminals
terminals <- setNames(as.numeric(uniprot_data$logFC), uniprot_data$external_gene_name)

# run PCSF with random noise

# time a pcsf run
start_time <- Sys.time()
subnet <- PCSF_rand(ppi, terminals, n = 50, r = 0.1, w = 2, b = 1, mu = 0.0005)
elapsed_time <- Sys.time() - start_time
print(elapsed_time)

plot.PCSF(subnet, node_label_cex = 15)







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


# load Fpocket data
af_drugability <- read.csv("../druggability_results/fpocket_druggability.csv")

# merge PCSF data with AF
PCSF_results <- merge(PCSF_master, af_drugability, by.x = "uniprot_gn_id", by.y = "uniprot_id")

# run citation scoring tool and read in
citation_scores <- read.csv("citation_scores.csv")

PCSF_results <- merge(PCSF_results, citation_scores, by.x = "external_gene_name", by.y = "gene_id")
PCSF_results <- na.omit(PCSF_results) # for some reason the merge created empty duplicates.
PCSF_results <- PCSF_results[order(-PCSF_results$druggability), ]
rownames(PCSF_results) <- NULL


# Jess wants to keep structure duplicates with highest druggability
filtered_results <- PCSF_results %>%
  group_by(external_gene_name) %>%
  filter(druggability == max(druggability))


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

