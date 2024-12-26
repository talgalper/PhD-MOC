library(randomForest)
library(caret)
library(progress)
library(pROC)


feature_matrix <- read.table("data/feature_matrix.txt", sep = "\t", header = T)

load("RData/TTD_data.RData")
target_info_df <- target_info_df[-1, ]

TTD_master <- data.frame(TARGETID = rownames(target_info_df),
                         GENENAME = target_info_df$GENENAME,
                         UNIPROID = target_info_df$UNIPROID,
                         TARGNAME = target_info_df$TARGNAME)

TTD_master <- merge(TTD_master, drug_info_df, by = "TARGETID", all = T)

TTD_master <- merge(TTD_master, drugDisease_data[, c(1:3)], by.x = "TTDDrugID", by.y = "TTDDRUID", all = T)
TTD_master <- unique(TTD_master)


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
# all_target_genes <- uniprot_ids_clean$uniprotswissprot
# all_target_genes <- na.omit(all_target_genes)

TTD_master <- merge(uniprot_ids_clean, TTD_master, by.x = "all_target_genes", by.y = "GENENAME", all = T)
TTD_master <- unique(TTD_master) # because i merged NAs

TTD_approved <- merge(unique(TTD_master[, c(1,2,4,10)]), TTD_approved, by = "TARGETID", all.y = T)
cancer_terms <- c("cancer", "carcinoma", "leukemia", "neoplasm", "metastases", "tumour")
pattern <- paste(cancer_terms, collapse = "|")
TTD_approved_cancer <- TTD_approved[grepl(pattern, TTD_approved$INDICATION, ignore.case = TRUE), ]
TTD_approved_cancer <- TTD_approved_cancer[!is.na(TTD_approved_cancer$INDICATION), ]

TTD_approved_cancer <- unique(TTD_approved_cancer$uniprotswissprot)
TTD_approved_cancer <- TTD_approved_cancer[TTD_approved_cancer != "" & !is.na(TTD_approved_cancer)]

feature_matrix$approved <- as.factor(ifelse(feature_matrix$Protein %in% TTD_approved_cancer, 1, 0))


TTD_clinical <- merge(unique(TTD_master[, c(1,2,4,10)]), TTD_clinical, by = "TARGETID", all.y = T)
cancer_terms <- c("cancer", "carcinoma", "leukemia", "neoplasm", "metastases")
pattern <- paste(cancer_terms, collapse = "|")
TTD_clinical_cancer <- TTD_clinical[grepl(pattern, TTD_clinical$INDICATION, ignore.case = TRUE), ]
TTD_clinical_cancer <- TTD_clinical_cancer[!is.na(TTD_clinical_cancer$INDICATION), ]

TTD_clinical_cancer <- unique(TTD_clinical_cancer$uniprotswissprot)
TTD_clinical_cancer <- TTD_clinical_cancer[TTD_clinical_cancer != "" & !is.na(TTD_clinical_cancer)]

feature_matrix$clinical <- as.factor(ifelse(feature_matrix$Protein %in% TTD_clinical_cancer, 1, 0))


# Split data into positive (approved drug targets) and negative (non-drug targets across any indication)
positive_set <- feature_matrix[feature_matrix$approved == 1, ]
negative_pool <- feature_matrix[!feature_matrix$Protein %in% all_target_genes, ]


# Initialize parameters
ntrees <- 1000  # Fixed number of trees
n_models <- 100  # Number of random forests for bagging
predictions <- matrix(0, nrow = nrow(feature_matrix), ncol = n_models)  # For storing predictions


pb <- progress_bar$new(format = "[:bar] :current/:total (:percent) eta: :eta", 
                       total = n_models)
# Train multiple RF models using random negative samples
for (i in 1:n_models) {
  set.seed(i)  # Different seed for each iteration
  # Resample negatives for each iteration
  negative_set <- negative_pool[sample(1:nrow(negative_pool), nrow(positive_set), replace = FALSE), ]
  training_data <- rbind(positive_set, negative_set)
  training_labels <- training_data$approved
  training_features <- training_data[, !names(training_data) %in% c("approved", "Protein", "clinical")]
  
  # Train Random Forest
  rf_model <- randomForest(
    x = training_features,
    y = as.factor(training_labels),
    mtry = round(sqrt(ncol(training_features))),
    ntree = ntrees,
    importance = TRUE)
  
  # Predict on the entire dataset
  predictions[, i] <- predict(rf_model, feature_matrix[, !names(feature_matrix) %in% c("approved", "Protein", "clinical")], type = "prob")[, 2]
  
  # save at the 100th and 1000th model 
  if (i == 100) {
    predictions_after_100 <- predictions
    predictions_after_100 <- predictions_after_100[1:100,1:100]
  }
  
  if (i == 1000) {
    predictions_after_1000 <- predictions
    predictions_after_1000 <- predictions_after_1000[1:1000,1:1000]
  }
  
  pb$tick()
}


# Average predictions across all models
final_predictions <- rowMeans(predictions)
# Attach predictions to the dataset
feature_matrix$Prediction_Score <- final_predictions

model_scores <- subset(feature_matrix, select = c("Protein", "Prediction_Score", "approved", "clinical"))

gene_names <- getBM(
  attributes = c("uniprotswissprot", "hgnc_symbol"),
  filters = "uniprotswissprot",
  values = model_scores$Protein,
  mart = ensembl)

model_scores <- merge(gene_names, model_scores, by.x = "uniprotswissprot", by.y = "Protein", all.y = T)
model_scores <- model_scores[order(-model_scores$Prediction_Score), ]


# Evaluate on validation set (if provided)
roc_curve <- roc(feature_matrix$clinical, feature_matrix$Prediction_Score)
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


save(roc_curve, importance_df, model_scores, feature_matrix, positive_set, negative_pool, file = "RData/updated_TTD_model.RData")
load("RData/updated_TTD_model.RData")


library(corrplot)

# Calculate correlations between features
# Exclude non-numeric columns like 'Protein', 'Label', or 'validation'
numeric_features <- feature_data[, sapply(feature_data, is.numeric)]
correlation_matrix <- cor(numeric_features, use = "pairwise.complete.obs")

# Plot the correlation matrix
corrplot(correlation_matrix, method = "color", type = "upper",
         tl.cex = 0.7, tl.col = "black", tl.srt = 45,
         title = "Feature Correlation Matrix", mar = c(0, 0, 1, 0))




