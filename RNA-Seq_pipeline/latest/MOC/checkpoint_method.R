library(tidyverse)

# checkpoint 1 is the top 100 differentially expressed genes
checkpoint_1 <- data.frame(up_reg = kylie_pcsf_master$external_gene_name[order(-kylie_pcsf_master$logFC)])
checkpoint_1 <- checkpoint_1[1:100, ]

# checkpoint 2 uses RRA to combine degree centrality and betweenness to get the top 100 genes "important" to the network
# genes common to this list and the checkpoint 1 list make up the checkpoint 2 file
rankings <- list(kylie_pcsf_master$external_gene_name[order(-kylie_pcsf_master$betweenness)],
                 kylie_pcsf_master$external_gene_name[order(-kylie_pcsf_master$degree_centrality)])
aggregate_ranks <- aggregateRanks(glist = rankings)
aggregate_ranks <- aggregate_ranks[1:100, ]

checkpoint_2 <- intersect(checkpoint_1, aggregate_ranks$Name)


# checkpoint 3 adds the druggability data from Fpocket and PocketMiner to be used to sort the checkpoint 2 list
drug_subset <- subset(kylie_pcsf_master, select = c("external_gene_name", "druggability", 
                                                    "num_drug_pockets", "max_hit", "num_hits"))

checkpoint_3 <- drug_subset[drug_subset$external_gene_name %in% checkpoint_2, ]


citation_scores <- read.csv("intermediate/MeSH_citation_scores.csv")
citation_subset <- citation_scores[citation_scores$gene_id %in% checkpoint_3$external_gene_name, ]

# citaiton scores are added in checkpoint 4
checkpoint_4 <- merge(checkpoint_3, citation_subset, by.x = "external_gene_name", by.y = "gene_id")
colnames(checkpoint_4)[5] <- "num_cryp_pockets"
