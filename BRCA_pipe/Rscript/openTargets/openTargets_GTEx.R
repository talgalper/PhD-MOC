library(geneSynonym)
library(biomaRt)
library(tidyverse)
library(reshape2)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")



# targets from API query (profile)
OpenTargets <- read_tsv("OpenTargets_data/breast_carcinoma_known_drugs.tsv")
#OpenTargets <- subset(OpenTargets, select = c("Target Approved Symbol", "Disease Name", "Drug Name", "Action Type", "Mechanism of Action","Drug Type"))

OpenTargets_filtered <- unique(OpenTargets) # drop duplicate rows
OpenTargets_filtered <- OpenTargets_filtered[OpenTargets_filtered$`Drug Type` == "Small molecule", ] # subset small molecule drugs
OpenTargets_filtered <- OpenTargets_filtered[OpenTargets_filtered$Status %in% c("Active, not recruiting", "Completed", "Recruiting"), ]

# either filter by unique drug name or gene name
OpenTargets_unique <- OpenTargets_filtered[!duplicated(OpenTargets_filtered$`Target Approved Symbol`), ] 
OpenTargets_unique <- OpenTargets_filtered[!duplicated(OpenTargets_filtered$`Drug Name`), ]
#OpenTargets_unique <- subset(OpenTargets_unique, select = c("Phase", "Disease Name", "Drug Name", "Target Approved Symbol",
#                                                            "Target Approved Name", "Mechanism of Action", "Action Type"))

OpenTargets_unique$NCT_ID <- sub(".*(NCT\\d+).*", "\\1", OpenTargets_unique$URL) # isolate NCT IDs


# read in data from NCT
NCT_summaries <- read.csv("OpenTargets_data/NCT_summaries.csv")
colnames(NCT_summaries) <- c("NCT_ID", "Breif Summary")

NCT_OpenTargets <- merge(OpenTargets_unique, NCT_summaries, by = "NCT_ID")
NCT_OpenTargets <- NCT_OpenTargets[!duplicated(NCT_OpenTargets$`Drug Name`), ]
NCT_OpenTargets <- subset(NCT_OpenTargets, select = c("NCT_ID", "Phase", "Status", "Disease Name", "Drug Name",
                                                      "Mechanism of Action", "Action Type", "Target Approved Name", 
                                                      "Target Approved Symbol", "Drug Type", "Breif Summary"))


# read in manually tagged NCT data
OpenTargets_NCT <- read.csv("OpenTargets_data/OpenTargets_NCT.csv", row.names = 1)
rownames(OpenTargets_NCT) <- NULL
table(OpenTargets_NCT$Cancer.Type)
table(OpenTargets_NCT$Cancer.Type[is.na(OpenTargets_NCT$Subtype)])

# existing annotations
OpenTargets_NCT$Subtype <- ifelse(OpenTargets_NCT$Disease.Name == "triple-negative breast cancer", "TNBC", NA)
OpenTargets_NCT$Subtype <- ifelse(
  OpenTargets_NCT$Disease.Name == "HER2 Positive Breast Carcinoma", 
  "Her2", 
  OpenTargets_NCT$Subtype
)

# luminal A
OpenTargets_NCT$Subtype <- ifelse(
  OpenTargets_NCT$Cancer.Type %in% c("ER+/Her2-", "ER+/Her2-/FGFR-", "ER+/PR+/Her2-", "HR+/Her2-"), 
  "lumA",
  OpenTargets_NCT$Subtype
)

# luminal B
OpenTargets_NCT$Subtype <- ifelse(
  grepl("ER+", OpenTargets_NCT$Cancer.Type, fixed = TRUE) & grepl("Her2+", OpenTargets_NCT$Cancer.Type, fixed = TRUE),
  "lumB",
  OpenTargets_NCT$Subtype)

# TNBC
OpenTargets_NCT$Subtype <- ifelse(
  grepl("TNBC", OpenTargets_NCT$Cancer.Type, fixed = TRUE),
  "TNBC",
  OpenTargets_NCT$Subtype
)

# luminal A/B
OpenTargets_NCT$Subtype <- ifelse(
  OpenTargets_NCT$Cancer.Type %in% c("ER+", "AR+/ER+", "ER+/PR+", "HR+"), 
  "lumA/B",
  OpenTargets_NCT$Subtype
)

# luminal A/TNBC
OpenTargets_NCT$Subtype <- ifelse(
  OpenTargets_NCT$Cancer.Type == "Her2-", 
  "lumA/TNBC",
  OpenTargets_NCT$Subtype
)

# luminal B/Her2
OpenTargets_NCT$Subtype <- ifelse(
  OpenTargets_NCT$Cancer.Type == "Her2+", 
  "lumB/Her2",
  OpenTargets_NCT$Subtype
)

# Her2/TNBC
OpenTargets_NCT$Subtype <- ifelse(
  OpenTargets_NCT$Cancer.Type == "HR-", 
  "Her2/TNBC",
  OpenTargets_NCT$Subtype
)

# misc
OpenTargets_NCT$Subtype <- ifelse(
  !is.na(OpenTargets_NCT$Cancer.Type) & is.na(OpenTargets_NCT$Subtype),
  "misc",
  OpenTargets_NCT$Subtype)


OpenTargets_NCT_filtered <- OpenTargets_NCT[OpenTargets_NCT$Impact.On.Cacner.Progression != "no", ]


lumA_OpenTargets <- OpenTargets_NCT_filtered[grepl("lumA", OpenTargets_NCT_filtered$Subtype), ]
lumB_OpenTargets <-  OpenTargets_NCT_filtered[OpenTargets_NCT_filtered$Subtype %in% c("lumB", "lumA/B", "lumB/Her2"), ]
Her2_OpenTargets <- OpenTargets_NCT_filtered[grepl("Her2", OpenTargets_NCT_filtered$Subtype), ]
basal_OpenTargets <- OpenTargets_NCT_filtered[grepl("TNBC", OpenTargets_NCT_filtered$Subtype), ]
misc_OpenTargets <- OpenTargets_NCT_filtered[OpenTargets_NCT_filtered$Subtype %in% c(NA, "misc"), ]



lumA_rank <- read.csv("intermediate/LumA/filterByExp/GTEx/final_gene_counts.csv")
lumA_rank <- rownames_to_column(lumA_rank)
colnames(lumA_rank)[1] <- "lumA_rank"
lumA_rank <- subset(lumA_rank, select = c("lumA_rank", "external_gene_name"))

lumB_rank <- read.csv("intermediate/LumB/filterByExp/GTEx/final_gene_counts.csv")
lumB_rank <- rownames_to_column(lumB_rank)
colnames(lumB_rank)[1] <- "lumB_rank"
lumB_rank <- subset(lumB_rank, select = c("lumB_rank", "external_gene_name"))

Her2_rank <- read.csv("intermediate/Her2/filterByExp/GTEx/final_gene_counts.csv")
Her2_rank <- rownames_to_column(Her2_rank)
colnames(Her2_rank)[1] <- "Her2_rank"
Her2_rank <- subset(Her2_rank, select = c("Her2_rank", "external_gene_name"))

basal_rank <- read.csv("intermediate/basal/filterByExp/GTEx/final_gene_counts.csv")
basal_rank <- rownames_to_column(basal_rank)
colnames(basal_rank)[1] <- "basal_rank"
basal_rank <- subset(basal_rank, select = c("basal_rank", "external_gene_name"))

ranks <- merge(lumA_rank, OpenTargets_NCT_filtered, by.x = "external_gene_name", by.y = "Target.Approved.Symbol", all = T)
ranks <- subset(ranks, select = c("external_gene_name", "lumA_rank", "Subtype", "Phase", "Drug.Name", "Mechanism.of.Action", "Impact.On.Cacner.Progression"))
ranks <- merge(lumB_rank, ranks, by = "external_gene_name", all = T)
ranks <- merge(Her2_rank, ranks, by = "external_gene_name", all = T)
ranks <- merge(basal_rank, ranks, by = "external_gene_name", all = T)

ranks <- ranks[!is.na(ranks$Drug.Name), ]

ranks <- merge(OpenTargets_NCT[c(1, 5)], ranks, by = "Drug.Name", all.y = T)

ranks <- ranks[!is.na(ranks$lumA_rank) | !is.na(ranks$lumB_rank) | !is.na(ranks$Her2_rank) | !is.na(ranks$basal_rank), ]


# proportion of drugs per subtype found






### merge with pipeline data
IDs <- getBM(attributes = c("external_gene_name", "ensembl_gene_id", "description", "uniprot_gn_id"), 
             filters = "external_gene_name", 
             values = OpenTargets_NCT_filtered$Target.Approved.Symbol, 
             mart = ensembl)
IDs$description <- gsub("\\s*\\[.*?\\]", "", IDs$description)

targets <- merge(IDs, OpenTargets_NCT_filtered, by.x = "external_gene_name", by.y = "Target.Approved.Symbol")



# add Fpocket druggability scores
af_drugability <- read.csv("../druggability_results/fpocket_druggability.csv")
targets <- merge(targets, af_drugability, by.x = "uniprot_gn_id", by.y = "uniprot_id", all.x = T)

# add PocketMiner scores
pocketminer_data <- read.csv("../pocketminer/results/pocketminer_results_3.0.csv")
targets <- merge(targets, pocketminer_data, by.x = "uniprot_gn_id", by.y = "uniprot_id", all.x = T)

load("RData/DE_results_master.RData")

# add logFC scores
lumA_hits <- DE_results$GTEx_lumA$hits
lumA_hits <- lumA_hits[lumA_hits$PValue <= 0.05, ]
lumA_hits <- subset(lumA_hits, select = c("gene_id", "logFC"))
colnames(lumA_hits)[2] <- "lumA_logFC"
targets <- merge(targets, lumA_hits, by.x = "ensembl_gene_id", by.y = "gene_id", all.x = T)

lumB_hits <- DE_results$GTEx_lumB[[1]]
lumB_hits <- lumB_hits[lumB_hits$PValue <= 0.05, ]
lumB_hits <- subset(lumB_hits, select = c("gene_id", "logFC"))
colnames(lumB_hits)[2] <- "lumB_logFC"
targets <- merge(targets, lumB_hits, by.x = "ensembl_gene_id", by.y = "gene_id", all.x = T)

Her2_hits <- DE_results$GTEx_Her2[[1]]
Her2_hits <- Her2_hits[Her2_hits$PValue <= 0.05, ]
Her2_hits <- subset(Her2_hits, select = c("gene_id", "logFC"))
colnames(Her2_hits)[2] <- "Her2_logFC"
targets <- merge(targets, Her2_hits, by.x = "ensembl_gene_id", by.y = "gene_id", all.x = T)

basal_hits <- DE_results$GTEx_basal[[1]]
basal_hits <- basal_hits[basal_hits$PValue <= 0.05, ]
basal_hits <- subset(basal_hits, select = c("gene_id", "logFC"))
colnames(basal_hits)[2] <- "basal_logFC"
targets <- merge(targets, basal_hits, by.x = "ensembl_gene_id", by.y = "gene_id", all.x = T)

# add centrality data
LumA_centrality <- read.csv("intermediate/LumA/filterByExp/GTEx/PCSF_output.csv")
colnames(LumA_centrality)[4] <- "lumA_centrality"
LumA_centrality <- subset(LumA_centrality, select = c("gene_id", "lumA_centrality"))
targets <- merge(targets, LumA_centrality, by.x = "external_gene_name", by.y = "gene_id", all.x = T)

lumB_centrality <- read.csv("intermediate/lumB/filterByExp/GTEx/PCSF_output.csv")
colnames(lumB_centrality)[4] <- "lumB_centrality"
lumB_centrality <- subset(lumB_centrality, select = c("gene_id", "lumB_centrality"))
targets <- merge(targets, lumB_centrality, by.x = "external_gene_name", by.y = "gene_id", all.x = T)

Her2_centrality <- read.csv("intermediate/Her2/filterByExp/GTEx/PCSF_output.csv")
colnames(Her2_centrality)[4] <- "Her2_centrality"
Her2_centrality <- subset(Her2_centrality, select = c("gene_id", "Her2_centrality"))
targets <- merge(targets, Her2_centrality, by.x = "external_gene_name", by.y = "gene_id", all.x = T)

basal_centrality <- read.csv("intermediate/basal/filterByExp/GTEx/PCSF_output.csv")
colnames(basal_centrality)[4] <- "basal_centrality"
basal_centrality <- subset(basal_centrality, select = c("gene_id", "basal_centrality"))
targets <- merge(targets, basal_centrality, by.x = "external_gene_name", by.y = "gene_id", all.x = T)

# add rank data
lumA_rank <- read.csv("intermediate/LumA/filterByExp/GTEx/final_gene_counts.csv")
lumA_rank <- rownames_to_column(lumA_rank)
colnames(lumA_rank)[1] <- "lumA_rank"
lumA_rank <- subset(lumA_rank, select = c("lumA_rank", "external_gene_name"))
targets <- merge(targets, lumA_rank, by.x = "external_gene_name", by.y = "external_gene_name", all.x = T)

lumB_rank <- read.csv("intermediate/LumB/filterByExp/GTEx/final_gene_counts.csv")
lumB_rank <- rownames_to_column(lumB_rank)
colnames(lumB_rank)[1] <- "lumB_rank"
lumB_rank <- subset(lumB_rank, select = c("lumB_rank", "external_gene_name"))
targets <- merge(targets, lumB_rank, by.x = "external_gene_name", by.y = "external_gene_name", all.x = T)

Her2_rank <- read.csv("intermediate/Her2/filterByExp/GTEx/final_gene_counts.csv")
Her2_rank <- rownames_to_column(Her2_rank)
colnames(Her2_rank)[1] <- "Her2_rank"
Her2_rank <- subset(Her2_rank, select = c("Her2_rank", "external_gene_name"))
targets <- merge(targets, Her2_rank, by.x = "external_gene_name", by.y = "external_gene_name", all.x = T)

basal_rank <- read.csv("intermediate/basal/filterByExp/GTEx/final_gene_counts.csv")
basal_rank <- rownames_to_column(basal_rank)
colnames(basal_rank)[1] <- "basal_rank"
basal_rank <- subset(basal_rank, select = c("basal_rank", "external_gene_name"))
targets <- merge(targets, basal_rank, by.x = "external_gene_name", by.y = "external_gene_name", all.x = T)


# tidy up
filtered_targets <- targets[!duplicated(targets$uniprot_gn_id) | !duplicated(targets$druggability, fromLast = TRUE) | !duplicated(targets$Drug.Name), ]
filtered_targets <- subset(filtered_targets, select = c("external_gene_name", "ensembl_gene_id", "description", 
                                                        "uniprot_gn_id", "pocket", "druggability", "struct_score","max_hit",
                                                        "Drug.Name", "Action.Type", "Mechanism.of.Action", "Drug.Type", "Disease.Name",
                                                        "lumA_logFC", "lumB_logFC", "Her2_logFC", "basal_logFC",
                                                        "lumA_centrality", "lumB_centrality", "Her2_centrality", "basal_centrality",
                                                        "lumA_rank", "lumB_rank", "Her2_rank", "basal_rank"))
filtered_targets <- filtered_targets[order(filtered_targets$external_gene_name), ]

# remove external_gene_name duplicates without a druggability score
duplicate_rows <- duplicated(filtered_targets$external_gene_name) | duplicated(filtered_targets$external_gene_name, fromLast = TRUE)
for (i in which(duplicate_rows)) {
  group <- filtered_targets[filtered_targets$external_gene_name == filtered_targets$external_gene_name[i], ]
  if (any(!is.na(group$druggability))) {
    filtered_targets <- filtered_targets[!(filtered_targets$external_gene_name == filtered_targets$external_gene_name[i] & is.na(filtered_targets$druggability)), ]
  }
}


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


ranks2 <- filtered_targets[!is.na(filtered_targets$lumA_rank) | !is.na(filtered_targets$lumB_rank) | !is.na(filtered_targets$Her2_rank) | !is.na(filtered_targets$basal_rank), ]
ranks2 <- ranks2[!duplicated(ranks2), ]

rownames(ranks2) <- NULL

ranks2 <- subset(ranks2, select = c("external_gene_name", "Drug.Name", "Mechanism.of.Action", "Action.Type",
                                    "lumA_rank", "lumB_rank", "Her2_rank", "basal_rank"))







### displays ranks across pipeline instead of subtypes ###
LumA_TCGA_norm <- read.csv("intermediate/LumA/filterByExp/GTEx/final_gene_counts.csv")
LumA_TCGA_norm <- rownames_to_column(LumA_TCGA_norm)
colnames(LumA_TCGA_norm)[1] <- "lumA_rank"
LumA_TCGA_norm <- subset(LumA_TCGA_norm, select = c("lumA_rank", "external_gene_name"))

LumB_TCGA_norm <- read.csv("intermediate/LumB/filterByExp/GTEx/final_gene_counts.csv")
LumB_TCGA_norm <- rownames_to_column(LumB_TCGA_norm)
colnames(LumB_TCGA_norm)[1] <- "lumB_rank"
LumB_TCGA_norm <- subset(LumB_TCGA_norm, select = c("lumB_rank", "external_gene_name"))

Her2_TCGA_norm <- read.csv("intermediate/Her2/filterByExp/GTEx/final_gene_counts.csv")
Her2_TCGA_norm <- rownames_to_column(Her2_TCGA_norm)
colnames(Her2_TCGA_norm)[1] <- "Her2_rank"
Her2_TCGA_norm <- subset(Her2_TCGA_norm, select = c("Her2_rank", "external_gene_name"))

basal_TCGA_norm <- read.csv("intermediate/basal/filterByExp/GTEx/final_gene_counts.csv")
basal_TCGA_norm <- rownames_to_column(basal_TCGA_norm)
colnames(basal_TCGA_norm)[1] <- "basal_rank"
basal_TCGA_norm <- subset(basal_TCGA_norm, select = c("basal_rank", "external_gene_name"))



LumA_GTEx <- read.csv("intermediate/LumA/filterByExp/GTEx/GTEx/final_gene_counts.csv")
LumA_GTEx <- rownames_to_column(LumA_GTEx)
colnames(LumA_GTEx)[1] <- "lumA_rank"
LumA_GTEx <- subset(LumA_GTEx, select = c("lumA_rank", "external_gene_name"))

LumB_GTEx <- read.csv("intermediate/LumB/filterByExp/GTEx/GTEx/final_gene_counts.csv")
LumB_GTEx <- rownames_to_column(LumB_GTEx)
colnames(LumB_GTEx)[1] <- "lumB_rank"
LumB_GTEx <- subset(LumB_GTEx, select = c("lumB_rank", "external_gene_name"))

Her2_GTEx <- read.csv("intermediate/Her2/filterByExp/GTEx/GTEx/final_gene_counts.csv")
Her2_GTEx <- rownames_to_column(Her2_GTEx)
colnames(Her2_GTEx)[1] <- "Her2_rank"
Her2_GTEx <- subset(Her2_GTEx, select = c("Her2_rank", "external_gene_name"))

basal_GTEx <- read.csv("intermediate/basal/filterByExp/GTEx/GTEx/final_gene_counts.csv")
basal_GTEx <- rownames_to_column(basal_GTEx)
colnames(basal_GTEx)[1] <- "basal_rank"
basal_GTEx <- subset(basal_GTEx, select = c("basal_rank", "external_gene_name"))



LumA_paired <- read.csv("intermediate/paired/LumA/final_gene_counts.csv")
LumA_paired <- rownames_to_column(LumA_paired)
colnames(LumA_paired)[1] <- "lumA_rank"
LumA_paired <- subset(LumA_paired, select = c("lumA_rank", "external_gene_name"))

LumB_paired <- read.csv("intermediate/paired/LumB/final_gene_counts.csv")
LumB_paired <- rownames_to_column(LumB_paired)
colnames(LumB_paired)[1] <- "lumB_rank"
LumB_paired <- subset(LumB_paired, select = c("lumB_rank", "external_gene_name"))

Her2_paired <- read.csv("intermediate/paired/Her2/final_gene_counts.csv")
Her2_paired <- rownames_to_column(Her2_paired)
colnames(Her2_paired)[1] <- "Her2_rank"
Her2_paired <- subset(Her2_paired, select = c("Her2_rank", "external_gene_name"))

basal_paired <- read.csv("intermediate/paired/basal/final_gene_counts.csv")
basal_paired <- rownames_to_column(basal_paired)
colnames(basal_paired)[1] <- "basal_rank"
basal_paired <- subset(basal_paired, select = c("basal_rank", "external_gene_name"))


lumA_common <- merge(LumA_TCGA_norm, LumA_GTEx, by = "external_gene_name", all = T)
lumA_common <- merge(lumA_common, LumA_paired, by = "external_gene_name", all = T)
colnames(lumA_common) <- c("external_gene_name", "TCGA_norm", "GTEx", "paired")
lumA_common[2:4] <- sapply(lumA_common[2:4], as.integer)

lumB_common <- merge(LumB_TCGA_norm, LumB_GTEx, by = "external_gene_name", all = T)
lumB_common <- merge(lumB_common, LumB_paired, by = "external_gene_name", all = T)
colnames(lumB_common) <- c("external_gene_name", "TCGA_norm", "GTEx", "paired")
lumB_common[2:4] <- sapply(lumB_common[2:4], as.integer)

Her2_common <- merge(Her2_TCGA_norm, Her2_GTEx, by = "external_gene_name", all = T)
Her2_common <- merge(Her2_common, Her2_paired, by = "external_gene_name", all = T)
colnames(Her2_common) <- c("external_gene_name", "TCGA_norm", "GTEx", "paired")
Her2_common[2:4] <- sapply(Her2_common[2:4], as.integer)

basal_common <- merge(basal_TCGA_norm, basal_GTEx, by = "external_gene_name", all = T)
basal_common <- merge(basal_common, basal_paired, by = "external_gene_name", all = T)
colnames(basal_common) <- c("external_gene_name", "TCGA_norm", "GTEx", "paired")
basal_common[2:4] <- sapply(basal_common[2:4], as.integer)


lumA_common <- merge(lumA_common, OpenTargets_NCT_filtered, by.x = "external_gene_name", by.y = "Target.Approved.Symbol")
lumA_common <- subset(lumA_common, select = c("external_gene_name", "TCGA_norm", "GTEx", "paired", 
                                              "Drug.Name", "Action.Type", "Subtype"))

lumB_common <- merge(lumB_common, OpenTargets_NCT_filtered, by.x = "external_gene_name", by.y = "Target.Approved.Symbol")
lumB_common <- subset(lumB_common, select = c("external_gene_name", "TCGA_norm", "GTEx", "paired", 
                                              "Drug.Name", "Action.Type", "Subtype"))

Her2_common <- merge(Her2_common, OpenTargets_NCT_filtered, by.x = "external_gene_name", by.y = "Target.Approved.Symbol")
Her2_common <- subset(Her2_common, select = c("external_gene_name", "TCGA_norm", "GTEx", "paired", 
                                              "Drug.Name", "Action.Type", "Subtype"))

basal_common <- merge(basal_common, OpenTargets_NCT_filtered, by.x = "external_gene_name", by.y = "Target.Approved.Symbol")
basal_common <- subset(basal_common, select = c("external_gene_name", "TCGA_norm", "GTEx", "paired", 
                                                "Drug.Name", "Action.Type", "Subtype"))











