

# run updated training and test data
library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
feature_matrix <- read.table("data/feature_matrix.txt", sep = "\t", header = T)
training_data <- data_sets_from_TTD(ensembl)

RF_results <- RF_bagging(training_data$feature_matrix, 
                         positive_set = training_data$positive_set, negative_pool = training_data$negative_pool, 
                         n_models = 1000, 
                         ntrees = 1000,
                         track_iterations = T, 
                         parallel = T, 
                         tuning = T, 
                         set_mtry_tuning_grid = seq(2, 10, by = 1))

