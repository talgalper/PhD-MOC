#### calculate overall druggability ####
targets <- read.csv("data_general/target_all_dbs.csv")
targets <- subset(targets, select = c("drugName", "drugBank_target", "Uniprot_ID", "ensembl_gene_id"))
targets <- targets[!duplicated(targets$Uniprot_ID), ]

Fpocket_scores <- read.csv("Fpocket/results_2024.05/fpocket_druggability.csv")
#Fpocket_scores <- merge(targets, Fpocket_scores, by.x = "Uniprot_ID", by.y = "uniprot_id", all.x = T)
#Fpocket_scores <- Fpocket_scores[order(Fpocket_scores$drugName), ]

PocketMiner_scores <- read.csv("PocketMiner/results/pocketminer_results_4.0.csv")
colnames(PocketMiner_scores)[3] <- "cryptic_pocket"


scores <- merge(targets, Fpocket_scores, by.x = "Uniprot_ID", by.y = "uniprot_id", all.x = T)
scores <- merge(scores, PocketMiner_scores, by.x = "Uniprot_ID", by.y = "uniprot_id", all.x = T)
scores <- subset(scores, select = c("Uniprot_ID", "drugName", "drugBank_target", "ensembl_gene_id", "druggability", "cryptic_pocket", "struct_score"))
scores <- scores[order(scores$drugName), ]

FDS <- scores$druggability
PCS <- scores$cryptic_pocket

# Set the gamma value
gamma <- 0.7
# Calculate the Adjusted Overall Druggability Score (ODS)
Adjusted_ODS <- gamma * pmax(FDS, PCS) + (1 - gamma) * pmin(FDS, PCS)

scores$ODS <- Adjusted_ODS




## minimum threshold
gamma <- 0.7

# Minimum threshold function
threshold_combination <- function(FDS, PCS, gamma, threshold = 0.8) {
  max_score <- pmax(FDS, PCS)
  min_score <- pmin(FDS, PCS)
  adjusted_score <- ifelse(max_score < threshold, gamma * max_score + (1 - gamma) * min_score, max_score)
  return(adjusted_score)
}

# Calculate the Adjusted Overall Druggability Score (ODS)
Adjusted_ODS <- threshold_combination(FDS, PCS, gamma)

# Display the results
scores$ODS <- Adjusted_ODS













DrugBank <- read.csv("DrugBank/DrugBank_targets_ENS.csv")
DrugBank <- merge(DrugBank, Fpocket_scores, by.x = "Uniprot_ID", by.y = "uniprot_id", all.x = T)
DrugBank <- merge(DrugBank, PocketMiner_scores, by.x = "Uniprot_ID", by.y = "uniprot_id", all.x = T)
DrugBank <- na.omit(DrugBank)

FDS <- DrugBank$druggability
PCS <- DrugBank$max_hit

gamma <- 0.8
Adjusted_ODS <- gamma * pmax(FDS, PCS) + (1 - gamma) * pmin(FDS, PCS)
DrugBank$ODS <- Adjusted_ODS
