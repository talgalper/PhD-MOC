### Make a data frame containing all information for known breast cancer subtypes ###


library(geneSynonym)
library(biomaRt)
library(tidyverse)
library(reshape2)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# get all synonyms for target terms
synon <- humanSyno(c("ESR1", "PGR", "ERBB2", "CDK4", "CDK6", 
                     "PI3K", "MTOR", "FGFR1", "FGFR2", "FGFR3", 
                     "FGFR4", "AKT", "ERK", "SRC", "PARP", 
                     "PD-L1", "MEK", "ERBB3", "AR", "TP53",
                     "MELK", "TOPK"))

# get terms from open targets platform
open_targets <- read.table("OT-MONDO_0007254-associated-targets-16_04_2024-v24_03.tsv", sep = "\t")
colnames(open_targets) <- open_targets[1, ]
open_targets <- open_targets[-1, ]
synon <- humanSyno(open_targets$symbol[1:50])


synon <- melt(synon)
colnames(synon) <- c("synonym", "NCBI_gene", "input_term")

# add addtional identifiers
IDs <- getBM(attributes = c("external_gene_name", "ensembl_gene_id", "description", "uniprot_gn_id"), 
                     filters = "external_gene_name", 
                     values = synon$synonym, 
                     mart = ensembl)
IDs$description <- gsub("\\s*\\[.*?\\]", "", IDs$description)
targets <- merge(synon, IDs,  by.x = "synonym", by.y = "external_gene_name")
targets <- subset(targets, select = c("input_term", "synonym", "ensembl_gene_id", 
                                      "description", "NCBI_gene", "uniprot_gn_id"))

# examine unrecognised synonyms
failed_terms <- merge(synon, IDs,  by.x = "synonym", by.y = "external_gene_name", all = T)

# add Fpocket druggability scores
af_drugability <- read.csv("../druggability_results/fpocket_druggability.csv")
targets <- merge(targets, af_drugability, by.x = "uniprot_gn_id", by.y = "uniprot_id", all.x = T)

# add PocketMiner scores
pocketminer_data <- read.csv("../pocketminer/results/pocketminer_results_3.0.csv")
targets <- merge(targets, pocketminer_data, by.x = "uniprot_gn_id", by.y = "uniprot_id", all.x = T)

# add logFC scores
lumA_hits <- read.csv("intermediate/paired/LumA/DE_results.csv")
lumA_hits <- lumA_hits[lumA_hits$PValue <= 0.05, ]
lumA_hits <- subset(lumA_hits, select = c("X", "logFC"))
colnames(lumA_hits)[2] <- "lumA_logFC"
targets <- merge(targets, lumA_hits, by.x = "ensembl_gene_id", by.y = "X", all.x = T)

lumB_hits <- read.csv("intermediate/paired/LumB/DE_results.csv")
lumB_hits <- lumB_hits[lumB_hits$PValue <= 0.05, ]
lumB_hits <- subset(lumB_hits, select = c("X", "logFC"))
colnames(lumB_hits)[2] <- "lumB_logFC"
targets <- merge(targets, lumB_hits, by.x = "ensembl_gene_id", by.y = "X", all.x = T)

Her2_hits <- read.csv("intermediate/paired/Her2/DE_results.csv")
Her2_hits <- Her2_hits[Her2_hits$PValue <= 0.05, ]
Her2_hits <- subset(Her2_hits, select = c("gene_id", "logFC"))
colnames(Her2_hits)[2] <- "Her2_logFC"
targets <- merge(targets, Her2_hits, by.x = "ensembl_gene_id", by.y = "gene_id", all.x = T)

basal_hits <- read.csv("intermediate/paired/basal/DE_results.csv")
basal_hits <- basal_hits[basal_hits$PValue <= 0.05, ]
basal_hits <- subset(basal_hits, select = c("gene_id", "logFC"))
colnames(basal_hits)[2] <- "basal_logFC"
targets <- merge(targets, basal_hits, by.x = "ensembl_gene_id", by.y = "gene_id", all.x = T)

# add centrality data
LumA_centrality <- read.csv("intermediate/paired/LumA/PCSF_output.csv")
colnames(LumA_centrality)[4] <- "lumA_centrality"
LumA_centrality <- subset(LumA_centrality, select = c("gene_id", "lumA_centrality"))
targets <- merge(targets, LumA_centrality, by.x = "synonym", by.y = "gene_id", all.x = T)

lumB_centrality <- read.csv("intermediate/paired/lumB/PCSF_output.csv")
colnames(lumB_centrality)[4] <- "lumB_centrality"
lumB_centrality <- subset(lumB_centrality, select = c("gene_id", "lumB_centrality"))
targets <- merge(targets, lumB_centrality, by.x = "synonym", by.y = "gene_id", all.x = T)

Her2_centrality <- read.csv("intermediate/paired/Her2/PCSF_output.csv")
colnames(Her2_centrality)[4] <- "Her2_centrality"
Her2_centrality <- subset(Her2_centrality, select = c("gene_id", "Her2_centrality"))
targets <- merge(targets, Her2_centrality, by.x = "synonym", by.y = "gene_id", all.x = T)

basal_centrality <- read.csv("intermediate/paired/basal/PCSF_output.csv")
colnames(basal_centrality)[4] <- "basal_centrality"
basal_centrality <- subset(basal_centrality, select = c("gene_id", "basal_centrality"))
targets <- merge(targets, basal_centrality, by.x = "synonym", by.y = "gene_id", all.x = T)

# add rank data
lumA_rank <- read.csv("intermediate/paired/LumA/final_gene_counts.csv")
lumA_rank <- rownames_to_column(lumA_rank)
colnames(lumA_rank)[1] <- "lumA_rank"
lumA_rank <- subset(lumA_rank, select = c("lumA_rank", "external_gene_name"))
targets <- merge(targets, lumA_rank, by.x = "synonym", by.y = "external_gene_name", all.x = T)

lumB_rank <- read.csv("intermediate/paired/LumB/final_gene_counts.csv")
lumB_rank <- rownames_to_column(lumB_rank)
colnames(lumB_rank)[1] <- "lumB_rank"
lumB_rank <- subset(lumB_rank, select = c("lumB_rank", "external_gene_name"))
targets <- merge(targets, lumB_rank, by.x = "synonym", by.y = "external_gene_name", all.x = T)

Her2_rank <- read.csv("intermediate/paired/Her2/final_gene_counts.csv")
Her2_rank <- rownames_to_column(Her2_rank)
colnames(Her2_rank)[1] <- "Her2_rank"
Her2_rank <- subset(Her2_rank, select = c("Her2_rank", "external_gene_name"))
targets <- merge(targets, Her2_rank, by.x = "synonym", by.y = "external_gene_name", all.x = T)

basal_rank <- read.csv("intermediate/paired/basal/final_gene_counts.csv")
basal_rank <- rownames_to_column(basal_rank)
colnames(basal_rank)[1] <- "basal_rank"
basal_rank <- subset(basal_rank, select = c("basal_rank", "external_gene_name"))
targets <- merge(targets, basal_rank, by.x = "synonym", by.y = "external_gene_name", all.x = T)


# tidy up
filtered_targets <- targets[!duplicated(targets$synonym) | !duplicated(targets$druggability, fromLast = TRUE), ]
filtered_targets <- subset(filtered_targets, select = c("input_term", "synonym", "ensembl_gene_id", "description", "NCBI_gene", 
                                      "uniprot_gn_id", "pocket", "druggability", "struct_score","max_hit",
                                      "lumA_logFC", "lumB_logFC", "Her2_logFC", "basal_logFC",
                                      "lumA_centrality", "lumB_centrality", "Her2_centrality", "basal_centrality",
                                      "lumA_rank", "lumB_rank", "Her2_rank", "basal_rank"))
filtered_targets <- filtered_targets[order(filtered_targets$input_term), ]

# remove synonym duplicates without a druggability score
# rerun this chunk twice for some reason
duplicate_rows <- duplicated(filtered_targets$synonym) | duplicated(filtered_targets$synonym, fromLast = TRUE)
for (i in which(duplicate_rows)) {
  group <- filtered_targets[filtered_targets$synonym == filtered_targets$synonym[i], ]
  if (any(!is.na(group$druggability))) {
    filtered_targets <- filtered_targets[!(filtered_targets$synonym == filtered_targets$synonym[i] & is.na(filtered_targets$druggability)), ]
  }
}


# alternate organisaiton
filtered_targets <- subset(filtered_targets, select = c("input_term", "synonym", "ensembl_gene_id", "description", "NCBI_gene", 
                                                        "uniprot_gn_id", "pocket", "druggability", "struct_score","max_hit",
                                                        "lumA_logFC", "lumA_centrality", "lumA_rank", "lumB_logFC", "lumB_centrality",
                                                        "lumB_rank", "Her2_logFC", "Her2_centrality", "Her2_rank", "basal_logFC",
                                                        "basal_centrality", "basal_rank"))


# flag shit logFC scores
filtered_targets$lumA_logFC <- ifelse(filtered_targets$lumA_logFC > -1 & filtered_targets$lumA_logFC < 1,
                                      paste(filtered_targets$lumA_logFC, "*"), 
                                      filtered_targets$lumA_logFC)
filtered_targets$lumB_logFC <- ifelse(filtered_targets$lumB_logFC > -1 & filtered_targets$lumB_logFC < 1,
                                      paste(filtered_targets$lumB_logFC, "*"), 
                                      filtered_targets$lumB_logFC)
filtered_targets$Her2_logFC <- ifelse(filtered_targets$Her2_logFC > -1 & filtered_targets$Her2_logFC < 1,
                                      paste(filtered_targets$Her2_logFC, "*"), 
                                      filtered_targets$Her2_logFC)
filtered_targets$basal_logFC <- ifelse(filtered_targets$basal_logFC > -1 & filtered_targets$basal_logFC < 1,
                                      paste(filtered_targets$basal_logFC, "*"), 
                                      filtered_targets$basal_logFC)

# flag shit druggability scores
filtered_targets$druggability <- ifelse(filtered_targets$druggability < 0.4,
                                      paste(filtered_targets$druggability, "*"), 
                                      filtered_targets$druggability)


# flag shit structure scores
filtered_targets$struct_score <- ifelse(filtered_targets$struct_score < 60,
                                        paste(filtered_targets$struct_score, "*"), 
                                        filtered_targets$struct_score)










