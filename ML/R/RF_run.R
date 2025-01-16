
# Run updated training and test data
library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
feature_matrix <- read.table("data/feature_matrix.txt", sep = "\t", header = T)
training_data <- data_sets_from_TTD(ensembl)

RF_results3 <- RF_bagging(training_data$feature_matrix, 
                         positive_set = training_data$positive_set, negative_pool = training_data$negative_pool, 
                         n_models = 1000, 
                         ntrees = 1000,
                         track_iterations = T, 
                         parallel = T, 
                         tuning = F, 
                         set_mtry_tuning_grid = seq(2, 10, by = 1))

save(RF_results, file = "results/RF_results_1000m_tunedParallel.RData")

# Run paper training data with tuning and random sampling
# Need to use updated negative pool because not provided by author
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









