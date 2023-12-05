library(tidyverse)
library(RobustRankAggreg)
library(biomaRt)
library(rDGIdb)


ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

kylie_pcsf_master <- read.csv("results/MOC_PCSF_master_unique.csv", row.names = 1)
pocketminer_data <- read.csv("../../../pocketminer/results/pocketminer_results_3.0.csv")
kylie_pcsf_master <- merge(kylie_pcsf_master, pocketminer_data, by = "ID")


# filter genes with logFC cutoff of 2
checkpoint_1 <- kylie_pcsf_master[kylie_pcsf_master$logFC >= 2 | kylie_pcsf_master$logFC <= -2, ]
checkpoint_1 <- subset(checkpoint_1, select = c("external_gene_name", "logFC", "betweenness", "degree_centrality", 
                                                "druggability", "num_drug_pockets", "max_hit", "num_hits"))
colnames(checkpoint_1)[7] <- "cryp_pocket"
colnames(checkpoint_1)[8] <- "num_cryp_pockets"


citation_scores <- read.csv("intermediate/citation_scores_3.0.csv")
citation_subset <- citation_scores[citation_scores$gene_id %in% checkpoint_1$external_gene_name, ]

# citaiton scores are added in checkpoint 4
checkpoint_2 <- merge(checkpoint_1, citation_subset, by.x = "external_gene_name", by.y = "gene_id")


# add gene description
description <- getBM(attributes = c("external_gene_name", "description"), 
                     filters = "external_gene_name", 
                     values = checkpoint_2$external_gene_name, 
                     mart = ensembl)

description$description <- gsub("\\s*\\[.*?\\]", "", description$description)


checkpoint_2 <- merge(description, checkpoint_2, by = "external_gene_name")



# network rank
network_rank <- list(kylie_pcsf_master$external_gene_name[order(-kylie_pcsf_master$betweenness)],
                     kylie_pcsf_master$external_gene_name[order(-kylie_pcsf_master$degree_centrality)])

network_RRA <- aggregateRanks(glist = network_rank)
network_RRA <- network_RRA[network_RRA$Name %in% checkpoint_2$external_gene_name, ]
colnames(network_RRA)[2] <- "network_rank"

# druggability rank
druggability_rank <- list(kylie_pcsf_master$external_gene_name[order(-kylie_pcsf_master$druggability)],
                          kylie_pcsf_master$external_gene_name[order(-kylie_pcsf_master$num_drug_pockets)])

druggability_RRA <- aggregateRanks(glist = druggability_rank)
druggability_RRA <- druggability_RRA[druggability_RRA$Name %in% checkpoint_2$external_gene_name, ]
colnames(druggability_RRA)[2] <- "druggability_rank"

# cryptic pocket rank
cryptic_rank <- list(kylie_pcsf_master$external_gene_name[order(-kylie_pcsf_master$max_hit)],
                          kylie_pcsf_master$external_gene_name[order(-kylie_pcsf_master$num_hits)])

cryptic_RRA <- aggregateRanks(glist = cryptic_rank)
cryptic_RRA <- cryptic_RRA[cryptic_RRA$Name %in% checkpoint_2$external_gene_name, ]
colnames(cryptic_RRA)[2] <- "cryptic_rank"

# merge together
checkpoint_3 <- merge(checkpoint_2, network_RRA, by.x = "external_gene_name", by.y = "Name")
checkpoint_3 <- merge(checkpoint_3, druggability_RRA, by.x = "external_gene_name", by.y = "Name")
checkpoint_3 <- merge(checkpoint_3, cryptic_RRA, by.x = "external_gene_name", by.y = "Name")


# drug database cross reference
DGIdb <- queryDGIdb(checkpoint_3$external_gene_name)
results <- byGene(DGIdb)



checkpoint_3 <- merge(checkpoint_3, subset(results, select = c("Gene", "DistinctDrugCount")), by.x = "external_gene_name", by.y = "Gene")


# subset a list of top hits
top_hits <- checkpoint_3[(checkpoint_3$druggability >= 0.7 | 
                           checkpoint_3$num_cryp_pockets >= 1) & 
                           checkpoint_3$logFC >= quantile(checkpoint_3$logFC, probs = 0.6)
                         , ]


top_hits[3:ncol(top_hits)] <- round(top_hits[3:ncol(checkpoint_3)], digits = 3)
detailed_results <- detailedResults(queryDGIdb(top_hits$external_gene_name))









### old version ###

library(tidyverse)
library(RobustRankAggreg)
library(biomaRt)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

kylie_pcsf_master <- read.csv("results/MOC_PCSF_master_unique.csv", row.names = 1)
pocketminer_data <- read.csv("../../../pocketminer/results/pocketminer_results_3.0.csv")
kylie_pcsf_master <- merge(kylie_pcsf_master, pocketminer_data, by = "ID")


# checkpoint 1 is the top 100 differentially expressed genes
#checkpoint_1 <- data.frame(up_reg = kylie_pcsf_master$external_gene_name[order(-kylie_pcsf_master$logFC)])
#checkpoint_1 <- checkpoint_1[1:100, ]

# checkpoint 1 alternative method with logFC cutoff of 2
checkpoint_1 <- kylie_pcsf_master[kylie_pcsf_master$logFC >= 2 | kylie_pcsf_master$logFC <= -2, ]


# checkpoint 2 uses RRA to combine degree centrality and betweenness to get the top 100 genes "important" to the network
# genes common to this list and the checkpoint 1 list make up the checkpoint 2 file
rankings <- list(kylie_pcsf_master$external_gene_name[order(-kylie_pcsf_master$betweenness)],
                 kylie_pcsf_master$external_gene_name[order(-kylie_pcsf_master$degree_centrality)])
aggregate_ranks <- aggregateRanks(glist = rankings)
aggregate_ranks <- aggregate_ranks[aggregate_ranks$Name %in% checkpoint_1$external_gene_name, ]
checkpoint_2 <- merge(checkpoint_1, aggregate_ranks, by.x = "external_gene_name", by.y = "Name")
checkpoint_2 <- checkpoint_2[order(checkpoint_2$Score), ]
checkpoint_2 <- checkpoint_2[, -ncol(checkpoint_2)]


checkpoint_3 <- subset(checkpoint_2, select = c("external_gene_name", "logFC", "druggability", 
                                                "num_drug_pockets", "max_hit", "num_hits"))





citation_scores <- read.csv("intermediate/citation_scores_3.0.csv")
citation_subset <- citation_scores[citation_scores$gene_id %in% checkpoint_3$external_gene_name, ]

# citaiton scores are added in checkpoint 4
checkpoint_4 <- merge(checkpoint_3, citation_subset, by.x = "external_gene_name", by.y = "gene_id")
colnames(checkpoint_4)[4] <- "cryp_pocket"
colnames(checkpoint_4)[5] <- "num_cryp_pockets"



description <- getBM(attributes = c("external_gene_name", "description"), 
                     filters = "external_gene_name", 
                     values = checkpoint_4$external_gene_name, 
                     mart = ensembl)

description$description <- gsub("\\s*\\[.*?\\]", "", description$description)


checkpoint_4 <- merge(description, checkpoint_4, by = "external_gene_name")



# subset top hit genes from checkpoint 4 based on custom parameters
top_hits <- checkpoint_4[checkpoint_4$druggability >= 0.7 & checkpoint_4$cryp_pocket >= 0.7, ]


