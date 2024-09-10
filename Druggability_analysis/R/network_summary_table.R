

pocket_scores <- merge(Fpocket_scores, PocketMiner_scores, by = c("ID", "uniprot_id"))
pocket_scores <- pocket_scores[order(-pocket_scores$druggability)]

# TP53, STAT3, MYC, TGFB1, TGFB2 and MET
non_druggable <- data.frame(uniprot_id = c("P04637", "P40763", "P01106", "P01137", "P61812", "P08581"), 
                            gene_id = c("TP53", "STAT3", "MYC", "TGFB1", "TGFB2", "MET"))
non_druggable <- merge(non_druggable, pocket_scores, by = "uniprot_id", all.x = T)
non_druggable <- subset(non_druggable, select = c("gene_id", "druggability", "cryptic_pocket"))



### create master DF with all network metrics for FDA targets ###


all_db_targets <- read.csv("data_general/target_all_dbs.csv")
FDA_drug_targets <- unique(all_db_targets$drugBank_target)


library(igraph)
library(tidyverse)

# load in and format HHnet subnet data
HHnet_result <- read_graph("../Hierarchical_HotNet/BRCA/results/hhnet_cluster1_netNeighs.graphml", format = "graphml")

HHnet_result <- data.frame(degree = degree(HHnet_result),
                 betweenness = betweenness(HHnet_result),
                 source = V(HHnet_result)$color)

HHnet_result$source <- ifelse(HHnet_result$source == "tomato", "subnet", "STRING")
HHnet_result <- rownames_to_column(HHnet_result)
HHnet_result <- HHnet_result[order(-HHnet_result$degree), ]
colnames(HHnet_result)[1] <- "gene_id"
rownames(HHnet_result) <- NULL


# load in a format PCSF subnet data
load("../BRCA_pipe/latest_run/RData/PCSF_subnet_FULL.RData")

clust <- components(subnet)
PCSF_result <- data.frame(gene_id = names(clust$membership), cluster = factor(clust$membership))
betweenness <- betweenness(subnet) 
centrality <- degree(subnet) 
PCSF_result$betweenness <- betweenness[as.character(PCSF_result$gene_id)]
PCSF_result$degree_centrality <- centrality[as.character(PCSF_result$gene_id)]
PCSF_result$betweenness <- as.integer(PCSF_result$betweenness)
PCSF_result$degree_centrality <- as.integer(PCSF_result$degree_centrality)
PCSF_result$prize <- V(subnet)$prize
PCSF_result$type <- V(subnet)$type

rownames(PCSF_result) <- 1:nrow(PCSF_result)

PCSF_result <- PCSF_result[order(-PCSF_result$degree_centrality), ]
rownames(PCSF_result) <- NULL



# load in and convert DE data to gene symbols
load("../WGCNA/BRCA/RData/DE_subset/dif_exp.RData")


# load in and format WGCNA data
load("../WGCNA/BRCA/RData/all_default/signed/all_kwithin.RData")
load("../WGCNA/BRCA/RData/all_default/signed/all_bwnet.RData")

load("../../../../OneDrive - RMIT University/PhD/large_git_files/WGCNA/TCGA_GTEx_filt_norm.RData") # Mac
load("../../../../Desktop/WGCNA_BRCA_large_files/TCGA_GTEx_filt_norm.RData") # Ubunut

load("../WGCNA/BRCA/RData/all_default/signed/venn_data.RData")

library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

WGCNAsigned_modules <- as.data.frame(bwnet$colors)
WGCNAsigned_modules <- rownames_to_column(WGCNAsigned_modules)
ensembl_converted <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description"), 
                           filters = "ensembl_gene_id", 
                           values = names(bwnet$colors), 
                           mart = ensembl)
ensembl_converted$description <- gsub("\\s*\\[.*?\\]", "", ensembl_converted$description)


unmapped <- ensembl_converted[ensembl_converted$external_gene_name == "", ]
unrecognised <- dif_exp[!dif_exp$gene_id %in% ensembl_converted$ensembl_gene_id, ]
ensembl_converted <- ensembl_converted[ensembl_converted$external_gene_name != "", ]
novel_transcripts <- unmapped[grep("novel transcript", unmapped$description), ]
novel_proteins <- unmapped[grep("novel protein", unmapped$description), ]
pseudogene <- unmapped[grep("pseudogene", unmapped$description), ]

WGCNAsigned_modules <- merge(WGCNAsigned_modules, ensembl_converted, by.x = "rowname", by.y = "ensembl_gene_id", all.x = T)
colnames(WGCNAsigned_modules) <- c("ensembl_id", "module", "external_gene_name", "description")

tumour_associated <- as.data.frame(tumour_associated)
tumour_associated <- merge(tumour_associated, WGCNAsigned_modules, by.x = "tumour_associated", by.y = "ensembl_id", all.x = T)


# combine all the data for the FDA drug targets
targets <- read.csv("data_general/target_all_dbs.csv")
targets <- unique(targets$drugBank_target)

master <- data.frame(target = targets)
master <- merge(master, WGCNAsigned_modules, by.x = "target", by.y = "external_gene_name", all.x = T)
master <- subset(master, select = c("target", "ensembl_id", "description", "module"))
master <- merge(master, HHnet_result, by.x = "target", by.y = "gene_id", all.x = T)
master <- master[, -c(6,7)]
colnames(master)[5] <- "HHnet_degree" 
master <- merge(master, PCSF_result, by.x = "target", by.y = "gene_id", all.x = T)
master <- master[, -c(6,7,9,10)]
colnames(master)[6] <- "PCSF_degree" 
master <- merge(master, dif_exp, by.x = "ensembl_id", by.y = "gene_id", all.x = T)
master <- master[, -c(8:11)]
master <- merge(master, kWithin, by.x = "ensembl_id", by.y = "row.names", all.x = T)
master <- master[, -c(8,10,11)]
master$tumour_associated <- ifelse(master$target %in% tumour_associated$external_gene_name, "yes", "no")

write.csv(master, "../../../../Desktop/temp.csv")
