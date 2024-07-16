library(tidyverse)
library(edgeR)
library(PCSF)

# Load data
load("RData/LumA/DE_data.RData")
load("RData/LumB/DE_data.RData")
load("RData/Her2/DE_data.RData")
load("RData/basal/DE_data.RData")
tumour_data <- cbind(LumA_unstranded, LumB_unstranded, Her2_unstranded, Basal_unstranded)

GTEx_data <- read.table("gene_reads_2017-06-05_v8_breast_mammary_tissue.gct", skip = 2)
colnames(GTEx_data) <- GTEx_data[1, ]
GTEx_data <- GTEx_data[-1, -1]
rownames(GTEx_data) <- NULL
GTEx_data$Name <- gsub("\\.\\d+", "", GTEx_data$Name)

# opt for having gene Ensembl IDs instead of gene names as rownames (same as TCGA)
GTEx_ENS <- column_to_rownames(GTEx_data, "Name")
GTEx_ENS <- GTEx_ENS[ , -1]
rownames <- rownames(GTEx_ENS)
GTEx_ENS <- as.data.frame(sapply(GTEx_ENS, as.numeric))
rownames(GTEx_ENS) <- rownames

rm(LumA_unstranded, LumB_unstranded, Her2_unstranded, Basal_unstranded, normal_unstranded)

data_filt <- filter_low_expr(disease_data = tumour_data,
                             control_data = GTEx_ENS)

DE_results <- DE_analysis(counts_matrix = data_filt$counts_filt,
                          sample_info = data_filt$sample_info)

DE_genes <- DE_results$dif_exp
hits <- DE_results$hits
topTags <- DE_results$toptags

save(DE_genes, hits, topTags, DE_gene_symbols, file = "RData/alt_pipe/DE_results.RData")

# read in drug data
OpenTargets <- read.csv("../BRCA_pipe/OpenTargets_data/OpenTargets_unique_drug.csv", row.names = 1)
OpenTargets_raw <- read_tsv("../BRCA_pipe/OpenTargets_data/breast_carcinoma_known_drugs.tsv")

NIH_targets <- read.table("../BRCA_pipe/NIH_BRCA_approved_drugs.txt", sep = "\t")
colnames(NIH_targets)[1] <- "approved_drugs"
NIH_targets$approved_drugs <- toupper(NIH_targets$approved_drugs)

approved_openTargets <- merge(NIH_targets, OpenTargets_raw, by.x = "approved_drugs", by.y = "Drug Name")
approved_openTargets <- approved_openTargets[!duplicated(approved_openTargets$approved_drugs) | !duplicated(approved_openTargets$`Target ID`), ]


temp <- approved_openTargets[approved_openTargets$`Target ID` %in% DE_genes$gene_id, ]
table(unique(approved_openTargets$`Target Approved Symbol`) %in% unique(temp$`Target Approved Symbol`))


write.table(DE_genes$gene_id, "../../../../Desktop/DE_genes.txt", row.names = F, col.names = F, quote = F)


### Try convert to gene IDS for STRING
### OUTCOME: Returned fewer results when input into STRING
#library(biomaRt)
#ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
#
#gene_id <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), 
#                 filters = "ensembl_gene_id", 
#                 values = DE_genes$gene_id, 
#                 mart = ensembl)
#
#gene_id_filt <- gene_id[gene_id$external_gene_name != "", ]
#gene_id_filt <- gene_id_filt[!duplicated(gene_id_filt$external_gene_name), ]
#
#
#missing_ids <- gene_id[gene_id$external_gene_name == "", ]
#unrecognised_genes <- DE_genes[!DE_genes$gene_id %in% gene_id$ensembl_gene_id, ]
#
#total_missing <- c(missing_ids$ensembl_gene_id, unrecognised_genes$gene_id)
#total_missing <- DE_genes[DE_genes$gene_id %in% total_missing, ]
#
#missing_genes_convert <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "uniprot_gn_id", "description"), 
#                                 filters = "ensembl_gene_id", 
#                                 values = total_missing$gene_id, 
#                                 mart = ensembl)
#
#missing_genes_convert <- merge(missing_genes_convert, total_missing, by.x = "ensembl_gene_id", by.y = "gene_id", all = T)
#missing_genes_convert <- missing_genes_convert[!duplicated(missing_genes_convert$ensembl_gene_id), ] # remove duplicate IDs casused by uniprot IDs
#
## number of unrecognised terms
#table(is.na(missing_genes_convert$uniprot_gn_id) & is.na(missing_genes_convert$description))
#
#novel_transcripts <- missing_genes_convert[grep("novel transcript", missing_genes_convert$description), ]
#novel_proteins <- missing_genes_convert[grep("novel protein", missing_genes_convert$description), ]
#pseudogene <- missing_genes_convert[grep("pseudogene", missing_genes_convert$description), ]
#
#
#write.table(gene_id_filt$external_gene_name, "../../../../Desktop/DE_genes.txt", row.names = F, col.names = F, quote = F)


string_edge_data <- read.table("intermediate/alt_pipe/STRING network (physical) default edge(+100).csv", header = T, sep = ",", stringsAsFactors = F)
ppi_list <- subset(string_edge_data, select = c("name", "stringdb..score"))
ppi_list <- ppi_list %>% 
  separate(name, sep = " ", into = c("node_1", "del", "node_2"))
ppi_list <- subset(ppi_list, select = c("node_1", "node_2", "stringdb..score"))
ppi_list$node_1 <- gsub(".*.\\.", "", ppi_list$node_1)
ppi_list$node_2 <- gsub(".*.\\.", "", ppi_list$node_2)

string_node_data <- read.table("intermediate/alt_pipe/STRING network (physical) default node(+100).csv", header = T, sep = ",", stringsAsFactors = F)
node_list <- subset(string_node_data, select = c("name", "display.name"))
node_list$name <- gsub(".*.\\.", "", node_list$name)
ppi_list$original_order <- seq_len(nrow(ppi_list))
merged_df <- merge(ppi_list, node_list, by.x = "node_1", by.y = "name", all.x = TRUE)
merged_df <- merge(merged_df, node_list, by.x = "node_2", by.y = "name", all.x = TRUE)
merged_df <- merged_df[order(merged_df$original_order), ]

final_df <- merged_df[, c("display.name.x", "display.name.y", "stringdb..score")]
colnames(final_df) <- c("node_1", "node_2", "score")

table(unique(approved_openTargets$`Target Approved Symbol`) %in% unique(string_node_data$display.name))


DE_gene_symbols <- merge(DE_genes, string_node_data, by.x = "gene_id", by.y = "query.term")
DE_gene_symbols <- subset(DE_gene_symbols, select = c("gene_id", "display.name", "logFC", "logCPM", 
                                                      "F", "PValue", "FDR"))


# set seed for reproducibility 
set.seed(1234)
# construct interactome
ppi <- construct_interactome(final_df)
# set terminals
terminals <- setNames(as.numeric(DE_gene_symbols$logFC), DE_gene_symbols$display.name)
# run PCSF with random noise
start_time <- Sys.time()
subnet <- PCSF_rand(ppi, terminals, n = 50, r = 0.1, w = 2, b = 1, mu = 0.0005)
elapsed_time <- Sys.time() - start_time
print(elapsed_time)

plot.PCSF(subnet, node_label_cex = 15)

save(subnet, file = "RData/alt_pipe/PCSF_subnet(+100).RData")

load("RData/alt_pipe/PCSF_subnet(+100).RData")

# extract cluster data
clust <- components(subnet)
df <- data.frame(gene_id = names(clust$membership), cluster = factor(clust$membership))
betweenness <- betweenness(subnet) 
centrality <- degree(subnet) 
df$betweenness <- betweenness[as.character(df$gene_id)]
df$degree_centrality <- centrality[as.character(df$gene_id)]
df$betweenness <- as.integer(df$betweenness)
df$degree_centrality <- as.integer(df$degree_centrality)
rownames(df) <- NULL

node_list <- subset(string_node_data, select = c("name", "query.term", "display.name", "stringdb..canonical.name"))
PCSF_df <- merge(node_list, df, by.x = "query.term", by.y = "gene_id", all.y = T)
PCSF_df <- PCSF_df[order(-PCSF_df$degree_centrality), ]

rm(ppi_list, string_edge_data, string_node_data, merged_df, node_list, df)

table(unique(approved_openTargets$`Target Approved Symbol`) %in% unique(PCSF_df_base$display.name))
table(unique(approved_openTargets$`Target Approved Symbol`) %in% unique(PCSF_df$display.name))


# add Fpocket data
af_drugability <- read.csv("../druggability_results/fpocket_druggability.csv")
PCSF_druggability <- merge(PCSF_df, af_drugability, by.x = "stringdb..canonical.name", by.y = "uniprot_id", all.x = T)
PCSF_druggability <- PCSF_druggability[,-8]
# keep structure duplicates with highest druggability
PCSF_druggability <- PCSF_druggability[order(-PCSF_druggability$druggability), ]
PCSF_druggability <- PCSF_druggability[!duplicated(PCSF_druggability$display.name), ]

# add cryptic pocket data
pocketminer_data <- read.csv("../pocketminer/results/pocketminer_results_3.0.csv")
PCSF_master <- merge(PCSF_druggability, pocketminer_data, by = "ID", all.x = T)
PCSF_master <- PCSF_master[order(-PCSF_master$degree_centrality), ]

# merge back with DE data
PCSF_master <- merge(PCSF_master, DE_genes)

# druggability score
FDS <- PCSF_master$druggability[!is.na(PCSF_master$druggability)]
PCS <- PCSF_master$max_hit[!is.na(PCSF_master$max_hit)]
gamma <- 0.7

Adjusted_ODS <- gamma * pmax(FDS, PCS) + (1 - gamma) * pmin(FDS, PCS)





### To do:
# read in new pocketminer data
# re-run PCSF with added STRING nodes (+100). fix query.term issue for final_df. 
#vneed to define node type as well steiner/terminal
# get PCSF iteration appearances. Could use attribute function from igraph.
# create df with all data
# adjust ranking until most of the drug targets come through
# determine which DE genes are missing from STRING & AF
# perform WGCNA to see if any of the missing DE genes are co-expressed/interact with known drug targets/drivers
# find a way to integrate the knowledge based network with the WGCNA results







