# re-creating results from back in the day: Rank sensitivity

library(data.table)
library(igraph)
library(progress)
library(tidyverse)


druggability <- read.csv("../Druggability_analysis/data_general/druggability_scores_annot.csv")
load("../ML/RData/full_fpocket_results.RData")
results_master$Score <- (results_master$Score - min(results_master$Score)) / (max(results_master$Score) - min(results_master$Score))
results_master <- results_master[order(-results_master$`Druggability Score`, -results_master$Score), ]
results_master <- results_master[!duplicated(results_master$uniprot_id), ]
druggability <- subset(druggability, select = c("external_gene_name", "druggability", "CP_score", "highest_score"))
druggability <- unique(druggability)

# PCSF
load("latest_run/RData/STN_filt/PCSF_results.RData")
PCSF <- df
rm(df)

# HHnet
HHnet <- read.csv("../Hierarchical_HotNet/BRCA/STN_filt/results/df_subnet.csv")
HHnet_enrich <- read.csv("../Hierarchical_HotNet/BRCA/STN_filt/results/df_subnetNeighs.csv")


# PubTator3 counts and get total counts per gene
PubTator <- fread("~/OneDrive - RMIT University/PhD/large_git_files/PubTator3/citaiton_counts_ognsmAnnot.csv") # mac
PubTator$tax_id <- as.character(PubTator$tax_id)
PubTator$entrezgene_id <- as.character(PubTator$entrezgene_id)
PubTator[, combined := ifelse(symbol == "", entrezgene_id, symbol)]
PubTator_counts <- PubTator[
  , .(counts = sum(count)),
  by = combined
]
PubTator_counts <- PubTator_counts[order(-PubTator_counts$counts), ]
colnames(PubTator_counts)[1] <- "symbol"

PCSF <- merge.data.table(as.data.table(PCSF), PubTator_counts, by.x = "external_gene_name", by.y = "symbol", all.x = T)
PCSF <- PCSF[order(-PCSF$degree_centrality), ]
PCSF <- merge(PCSF, druggability, by = "external_gene_name", all.x = T)
PCSF <- PCSF[complete.cases(PCSF), ]

HHnet <- merge.data.table(as.data.table(HHnet), PubTator_counts, by.x = "external_gene_name", by.y = "symbol", all.x = T)
HHnet <- HHnet[order(-HHnet$degree), ]
HHnet <- merge(HHnet, druggability, by = "external_gene_name", all.x = T)
HHnet <- HHnet[complete.cases(HHnet), ]
HHnet_enrich <- merge.data.table(as.data.table(HHnet_enrich), PubTator_counts, by.x = "external_gene_name", by.y = "symbol", all.x = T)
HHnet_enrich <- HHnet_enrich[order(-HHnet_enrich$degree), ]
HHnet_enrich <- merge(HHnet_enrich, druggability, by = "external_gene_name", all.x = T)
HHnet_enrich <- HHnet_enrich[complete.cases(HHnet_enrich), ]

# normalise scores
PCSF$counts_norm <- log10(PCSF$counts)
PCSF[, c(6:9,16)] <- lapply(PCSF[, c(6:9,16)], function(x) {
  (x - min(x)) / (max(x) - min(x))
})

# choose weights
betweeness_w <- 1/6
centrality_w <- 1/6
closeness_w <- 1/6
eigen_centrality_w <- 1/6
druggability_w <- 1/6
citation_w <- 1/6

# or

betweeness_w <- 0.15
centrality_w <- 0.25
closeness_w <- 0.10
eigen_centrality_w <- 0.15
druggability_w <- 0.25
citation_w <- 0.10

# Combine scores using custom weights
combined_score <- betweeness_w * PCSF$betweenness +
  centrality_w * PCSF$degree_centrality +
  closeness_w * PCSF$closeness +
  eigen_centrality_w * PCSF$eigen_centrality +
  druggability_w * PCSF$highest_score -
  citation_w * PCSF$counts_norm

PCSF$rank_score <- combined_score
PCSF <- PCSF[order(-PCSF$rank_score), ]
PCSF <- rownames_to_column(PCSF)


targets <- read.csv("../Druggability_analysis/data_general/target_all_dbs.csv")
targets <- targets[, c(2,4)]
targets <- targets[!duplicated(targets$ensembl_gene_id), ]


temp <- PCSF[PCSF$external_gene_name %in% targets$drugBank_target, ]


# repeat for HHnet
HHnet$counts_norm <- log10(HHnet$counts)
HHnet[, c(5:8,15)] <- lapply(HHnet[, c(5:8,15)], function(x) {
  (x - min(x)) / (max(x) - min(x))
})



# Combine scores using custom weights
combined_score <- betweeness_w * HHnet$betweenness +
  centrality_w * HHnet$degree +
  closeness_w * HHnet$closeness +
  eigen_centrality_w * HHnet$eigen_centrality +
  druggability_w * HHnet$highest_score -
  citation_w * HHnet$counts_norm

HHnet$rank_score <- combined_score
HHnet <- HHnet[order(-HHnet$rank_score), ]
HHnet <- rownames_to_column(HHnet)




# repeat for HHnet neighbour enriched
HHnet_enrich$counts_norm <- log10(HHnet_enrich$counts)
HHnet_enrich[, c(5:8,15)] <- lapply(HHnet_enrich[, c(5:8,15)], function(x) {
  (x - min(x)) / (max(x) - min(x))
})


# Combine scores using custom weights
combined_score <- betweeness_w * HHnet_enrich$betweenness +
  centrality_w * HHnet_enrich$degree +
  closeness_w * HHnet_enrich$closeness +
  eigen_centrality_w * HHnet_enrich$eigen_centrality +
  druggability_w * HHnet_enrich$highest_score -
  citation_w * HHnet_enrich$counts_norm

HHnet_enrich$rank_score <- combined_score
HHnet_enrich <- HHnet_enrich[order(-HHnet_enrich$rank_score), ]
HHnet_enrich <- rownames_to_column(HHnet_enrich)


targets <- read.csv("../Druggability_analysis/data_general/target_all_dbs.csv")
targets <- targets[, c(2,4)]
targets <- targets[!duplicated(targets$ensembl_gene_id), ]


temp <- HHnet_enrich[HHnet_enrich$external_gene_name %in% targets$drugBank_target, ]









