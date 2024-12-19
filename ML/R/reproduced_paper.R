library(randomForest)
library(caret)
library(readxl)
library(progress)
library(pROC)

feature_data <- read.table("data/feature_matrix.txt", sep = "\t", header = T)
approved_targets <- read_xls("data/approved_targets.xls", sheet = 1)
clinical_targets <- read_xls("data/approved_targets.xls", sheet = 2)


# get all uniprot IDs for all TTD target proteins
all_target_genes <- unique(TTD_master$GENENAME)
all_target_genes <- na.omit(all_target_genes)
all_target_genes <- strsplit(all_target_genes, ";")
all_target_genes <- unlist(all_target_genes)
all_target_genes <- unique(all_target_genes)
all_target_genes <- trimws(all_target_genes)
all_target_genes <- all_target_genes[all_target_genes != ""]

library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

uniprot_ids <- getBM(
  attributes = c("hgnc_symbol", "uniprotswissprot"),
  filters = "hgnc_symbol",
  values = all_target_genes,
  mart = ensembl)

# remove empty duplicates
duplicated_values <- names(table(uniprot_ids$hgnc_symbol)[table(uniprot_ids$hgnc_symbol) > 1])
uniprot_ids_clean <- uniprot_ids[!(uniprot_ids$hgnc_symbol %in% duplicated_values & uniprot_ids$uniprotswissprot == ""), ]

uniprot_ids_clean <- merge(as.data.frame(all_target_genes), uniprot_ids_clean, by.x = "all_target_genes", by.y = "hgnc_symbol", all.x = T)
all_target_genes <- uniprot_ids_clean$uniprotswissprot
all_target_genes <- na.omit(all_target_genes)



feature_data$Label <- as.factor(ifelse(feature_data$Protein %in% approved_targets$Protein, 1, 0))
feature_data$validation <- as.factor(ifelse(feature_data$Protein %in% clinical_targets$Protein, 1, 0))

# Split data into positive (approved drug targets) and negative (non-drug targets)
positive_set <- feature_data[feature_data$Label == 1, ]
negative_pool <- feature_data[feature_data$Label == 0 & !feature_data$Protein %in% all_target_genes, ]


# Initialize parameters
ntrees <- 1000  # Fixed number of trees
n_models <- 100  # Number of random forests for bagging
predictions <- matrix(0, nrow = nrow(feature_data), ncol = n_models)  # For storing predictions


pb <- progress_bar$new(format = "[:bar] :current/:total (:percent) eta: :eta", 
                       total = n_models)
# Train multiple RF models using random negative samples
for (i in 1:n_models) {
  set.seed(i)  # Different seed for each iteration
  # Resample negatives for each iteration
  negative_set <- negative_pool[sample(1:nrow(negative_pool), nrow(positive_set), replace = FALSE), ]
  training_data <- rbind(positive_set, negative_set)
  training_labels <- training_data$Label
  training_features <- training_data[, !names(training_data) %in% c("Label", "Protein", "validation")]
  
  # Train Random Forest
  rf_model <- randomForest(
    x = training_features,
    y = as.factor(training_labels),
    mtry = round(sqrt(ncol(training_features))),
    ntree = ntrees,
    importance = TRUE)
  
  # Predict on the entire dataset
  predictions[, i] <- predict(rf_model, feature_data[, !names(feature_data) %in% c("Label", "Protein", "validation")], type = "prob")[, 2]
  
  pb$tick()
}


# Average predictions across all models
final_predictions <- rowMeans(predictions)
# Attach predictions to the dataset
feature_data$Prediction_Score <- final_predictions

model_scores <- subset(feature_data, select = c("Protein", "Prediction_Score", "Label", "validation"))

gene_names <- getBM(
  attributes = c("uniprotswissprot", "hgnc_symbol"),
  filters = "uniprotswissprot",
  values = model_scores$Protein,
  mart = ensembl)

model_scores <- merge(gene_names, model_scores, by.x = "uniprotswissprot", by.y = "Protein", all.y = T)
model_scores <- model_scores[order(-model_scores$Prediction_Score), ]


# Evaluate on validation set (if provided)
# validation_set <- feature_data[feature_data$validation == 1, ]  # Assuming 'Validation' column marks validation samples
roc_curve <- roc(feature_data$validation, feature_data$Prediction_Score)
auc(roc_curve)

# Plot ROC Curve
plot(roc_curve, main = "ROC Curve for Random Forest Predictions")


# Extract feature importance from one RF model as an example
importance_scores <- importance(rf_model)
importance_df <- data.frame(Feature = rownames(importance_scores), Importance = importance_scores[, "MeanDecreaseGini"])
importance_df <- importance_df[order(importance_df$Importance, decreasing = TRUE), ]


ggplot(importance_df[1:(0.5*nrow(importance_df)), ], aes(x = reorder(Feature, -Importance), y = Importance)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Feature Importance",
    x = "Features",
    y = "Mean Decrease in Gini"
  ) +
  theme_minimal()


save(roc_curve, importance_df, model_scores, feature_data, positive_set, negative_pool, file = "RData/reproduced_paper.RData")


