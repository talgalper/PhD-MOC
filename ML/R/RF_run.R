
#####################  Run updated training and test data ##################### 
library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
feature_matrix <- read.table("data/feature_matrix.txt", sep = "\t", header = T)
training_data <- data_sets_from_TTD(feature_matrix, ensembl)

RF_results <- RF_bagging(training_data$feature_matrix, 
                         positive_set = training_data$positive_set, negative_pool = training_data$negative_pool, 
                         n_models = 1000, 
                         ntrees = 1000,
                         track_iterations = T, 
                         parallel = T, 
                         tuning = T, 
                         set_mtry_tuning_grid = seq(2, 10, by = 1))

save(RF_results, file = "results/RF_results_1000m_tunedParallel.RData")

##################### Run paper training data with tuning and random sampling ##################### 
##################### Need to use updated negative pool because not provided by author ############
library(readxl)
feature_data <- read.table("data/feature_matrix.txt", sep = "\t", header = T)
approved_targets <- read_xls("data/approved_targets.xls", sheet = 1)
clinical_targets <- read_xls("data/approved_targets.xls", sheet = 2)

feature_data$approved <- as.factor(ifelse(feature_data$Protein %in% approved_targets$Protein, 1, 0))
feature_data$clinical <- as.factor(ifelse(feature_data$Protein %in% clinical_targets$Protein, 1, 0))
positive_set <- feature_data[feature_data$approved == 1, ]

RF_results2 <- RF_bagging(feature_matrix = feature_data, 
                         positive_set = positive_set, negative_pool = training_data$negative_pool, 
                         n_models = 1000, 
                         ntrees = 1000,
                         track_iterations = T, 
                         parallel = T, 
                         tuning = T, 
                         set_mtry_tuning_grid = seq(2, 10, by = 1))

save(RF_results2, file = "results/RF_results_1000m_tunedParallel_reproPaper.RData")


#####################  Run with extended feature matrix ##################### 
feature_matrix <- read.table("data/feature_matrix.txt", sep = "\t", header = T)
druggability <- read.csv("../Druggability_analysis/data_general/druggability_scores_annot.csv")
load("RData/full_fpocket_results.RData")
load("RData/human_string_PPI_metrics.RData")

results_master <- results_master[,-c(22,23)]
results_master$Score <- (results_master$Score - min(results_master$Score)) / (max(results_master$Score) - min(results_master$Score))
results_master <- results_master[order(-results_master$`Druggability Score`, -results_master$Score), ]
results_master <- results_master[!duplicated(results_master$uniprot_id), ]

druggability <- subset(druggability, select = c("uniprot_gn_id", "CP_score", "highest_score"))
druggability <- unique(druggability)

#new_feature_matrix <- feature_matrix[, -c((ncol(feature_matrix)-5):ncol(feature_matrix))] # remove old centrality data
new_feature_matrix <- feature_matrix
new_feature_matrix <- merge(new_feature_matrix, results_master[, -c(1,23)], by.x = "Protein", by.y = "uniprot_id", all.x = T)
new_feature_matrix <- merge(new_feature_matrix, subset(druggability, select = c("uniprot_gn_id", "CP_score", "highest_score")), by.x = "Protein", by.y = "uniprot_gn_id", all.x = T)
new_feature_matrix <- merge(new_feature_matrix, string_df, by.x = "Protein", by.y = "ENSG", all.x = T)

# convert all missing values to 0
new_feature_matrix[is.na(new_feature_matrix)] <- 0

# Normalise only the specified columns
new_feature_matrix[c(93:97)] <- lapply(new_feature_matrix[c(93:97)], function(x) {
  (x - min(x)) / (max(x) - min(x))
})

training_data <- data_sets_from_TTD(new_feature_matrix, ensembl)

RF_results3 <- RF_bagging(training_data$feature_matrix, 
                         positive_set = training_data$positive_set, negative_pool = training_data$negative_pool, 
                         n_models = 1000, 
                         ntrees = 1000,
                         track_iterations = T, 
                         parallel = T, 
                         tuning = T, 
                         set_mtry_tuning_grid = seq(2, 10, by = 1))

