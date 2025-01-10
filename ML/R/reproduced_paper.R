library(randomForest)
library(caret)
library(readxl)
library(progress)
library(pROC)

feature_data <- read.table("data/feature_matrix.txt", sep = "\t", header = T)
approved_targets <- read_xls("data/approved_targets.xls", sheet = 1)
clinical_targets <- read_xls("data/approved_targets.xls", sheet = 2)


feature_data$Label <- as.factor(ifelse(feature_data$Protein %in% approved_targets$Protein, 1, 0))
feature_data$validation <- as.factor(ifelse(feature_data$Protein %in% clinical_targets$Protein, 1, 0))

# Split data into positive (approved drug targets) and negative (non-drug targets)
# need to load in updated set of negative targets
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

positive_set <- feature_data[feature_data$Label == 1, ]
negative_pool <- feature_data[feature_data$Label == 0 & !feature_data$Protein %in% all_target_genes, ]


# Initialize parameters
ntrees <- 1000  # Fixed number of trees
n_models <- 10000  # Number of random forests for bagging
predictions <- matrix(0, nrow = nrow(feature_data), ncol = n_models)  # For storing predictions
# for storing importance scores. remove protein and label cols
n_features <- ncol(feature_data) - 3
importance_scores <- matrix(0, nrow = n_features, ncol = n_models) 

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
  
  # Extract feature importance 
  importance <- importance(rf_model)
  importance_scores[, i] <- importance[, "MeanDecreaseGini"]
  
  pb$tick()
}

save(predictions, importance_scores, file = "~/Desktop/reproduced_paper_predictions.RData")
load("~/Desktop/reproduced_paper_predictions.RData")

# Average predictions across all models
final_predictions <- rowMeans(predictions)

final_predictions_after_100 <- predictions[, 1:100]
final_predictions_after_100 <- rowMeans(final_predictions_after_100)

final_predictions_after_1000 <- predictions[, 1:1000]
final_predictions_after_1000 <- rowMeans(final_predictions_after_1000)

# Attach predictions to the dataset
feature_data$Prediction_Score_100m <- final_predictions_after_100
feature_data$Prediction_Score_1000m <- final_predictions_after_1000
feature_data$Prediction_Score_10000m <- final_predictions

write.csv(feature_data, "results/reproduced_paper_featureMatrix.csv")
feature_data <- read.csv("results/reproduced_paper_featureMatrix.csv", row.names = 1)

# average feature importance across all models
importance_df <- data.frame(Features = colnames(feature_data)[!colnames(feature_data) %in% c("Label", "Protein", "validation", "Prediction_Score_100m", "Prediction_Score_1000m", "Prediction_Score_10000m")],
                            Importance_100m = rowMeans(importance_scores[, 1:100]),
                            Importance_1000m = rowMeans(importance_scores[, 1:1000]),
                            Importance_10000m = rowMeans(importance_scores[, 1:10000]))

write.csv(importance_df, "results/feature_importance/reproduced_paper_FI.csv")
importance_df <- read.csv("results/feature_importance/reproduced_paper_FI.csv", row.names = 1)


model_scores <- subset(feature_data, select = c("Protein", "Prediction_Score_100m", "Prediction_Score_1000m", "Prediction_Score_10000m", "Label", "validation"))

library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

gene_names <- getBM(
  attributes = c("uniprotswissprot", "hgnc_symbol"),
  filters = "uniprotswissprot",
  values = model_scores$Protein,
  mart = ensembl)

model_scores <- merge(gene_names, model_scores, by.x = "uniprotswissprot", by.y = "Protein", all.y = T)
model_scores <- model_scores[order(-model_scores$Prediction_Score_10000m), ]
rownames(model_scores) <- NULL

write.csv(model_scores, "results/hgnc_symbol/reproduced_paper_modelScores.csv")

# Evaluate on validation set (if provided)
# validation_set <- feature_data[feature_data$validation == 1, ]  
roc_curve_100 <- roc(feature_data$validation, feature_data$Prediction_Score_100m)
roc_curve_1000 <- roc(feature_data$validation, feature_data$Prediction_Score_1000m)
roc_curve_10000 <- roc(feature_data$validation, feature_data$Prediction_Score_10000m)

auc(roc_curve_100)
auc(roc_curve_1000)
auc(roc_curve_10000)


# Plot ROC Curve
plot(roc_curve, main = "ROC Curve for Random Forest Predictions")

# plot importance df
ggplot(importance_df[1:(0.5*nrow(importance_df)), ], aes(x = reorder(Feature, -Importance), y = Importance)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Feature Importance",
    x = "Features",
    y = "Mean Decrease in Gini"
  ) +
  theme_minimal()


save(roc_curve_100, roc_curve_1000, roc_curve_10000, importance_df, model_scores, feature_data, positive_set, negative_pool, file = "RData/reproduced_paper.RData")




# see if updated list of validated targets were predicted
load("RData/reproduced_paper.RData")
load("RData/TTD_data.RData")

target_info_df <- target_info_df[-1, ]

TTD_master <- data.frame(TARGETID = rownames(target_info_df),
                         GENENAME = target_info_df$GENENAME,
                         UNIPROID = target_info_df$UNIPROID,
                         TARGNAME = target_info_df$TARGNAME)

TTD_master <- merge(TTD_master, drug_info_df, by = "TARGETID", all = T)

TTD_master <- merge(TTD_master, drugDisease_data[, c(1:3)], by.x = "TTDDrugID", by.y = "TTDDRUID", all = T)
TTD_master <- unique(TTD_master)

TTD_master <- merge(TTD_master, drug_rawinfo_df[, c(1,5)], by.x = "TTDDrugID", by.y = "DRUG__ID", all = T)
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

TTD_approved <- merge(unique(TTD_master[, c(1,2,4,10,11)]), TTD_approved, by = "TARGETID", all.y = T)
cancer_terms <- c("cancer", "carcinoma", "leukemia", "neoplasm", "metastases", "tumour")
pattern <- paste(cancer_terms, collapse = "|")
TTD_approved_cancer <- TTD_approved[grepl(pattern, TTD_approved$INDICATION, ignore.case = TRUE), ]
TTD_approved_cancer <- TTD_approved_cancer[!is.na(TTD_approved_cancer$INDICATION), ]
TTD_approved_cancer <- TTD_approved_cancer[grepl(c("small molecul"), TTD_approved$DRUGTYPE, ignore.case = TRUE), ] # molecule spelt wrong on purpose
TTD_approved_cancer <- TTD_approved_cancer[!is.na(TTD_approved_cancer$DRUGTYPE), ]

TTD_approved_cancer <- unique(TTD_approved_cancer$uniprotswissprot)
TTD_approved_cancer <- TTD_approved_cancer[TTD_approved_cancer != "" & !is.na(TTD_approved_cancer)]

feature_data$updated_approved <- as.factor(ifelse(feature_data$Protein %in% TTD_approved_cancer, 1, 0))


TTD_clinical <- merge(unique(TTD_master[, c(1,2,4,10,11)]), TTD_clinical, by = "TARGETID", all.y = T)
cancer_terms <- c("cancer", "carcinoma", "leukemia", "neoplasm", "metastases")
pattern <- paste(cancer_terms, collapse = "|")
TTD_clinical_cancer <- TTD_clinical[grepl(pattern, TTD_clinical$INDICATION, ignore.case = TRUE), ]
TTD_clinical_cancer <- TTD_clinical_cancer[!is.na(TTD_clinical_cancer$INDICATION), ]
# looks like they did not use small molecule clinical trial drugs for validation
# TTD_clinical_cancer <- TTD_clinical_cancer[grepl(c("small molecul"), TTD_clinical_cancer$DRUGTYPE, ignore.case = TRUE), ] # molecule spelt wrong on purpose
# TTD_clinical_cancer <- TTD_clinical_cancer[!is.na(TTD_clinical_cancer$DRUGTYPE), ]

TTD_clinical_cancer <- unique(TTD_clinical_cancer$uniprotswissprot)
TTD_clinical_cancer <- TTD_clinical_cancer[TTD_clinical_cancer != "" & !is.na(TTD_clinical_cancer)]

feature_data$updated_clinical <- as.factor(ifelse(feature_data$Protein %in% TTD_clinical_cancer, 1, 0))



plot_data <- subset(feature_data, select = c("Protein", "Prediction_Score_100m", "Prediction_Score_1000m", "Prediction_Score_10000m", "Label", "updated_approved", "validation", "updated_clinical"))
plot_data <- plot_data[order(-plot_data$Prediction_Score_10000m), ]

plot_data$Group <- with(plot_data, ifelse(
  Label == 1, "Approved Targets",
  ifelse(updated_approved == 1, "Updated Approved",
         ifelse(validation == 1, "Clinical Targets",
                ifelse(updated_clinical == 1, "Updated Clinical", "Not Targets")
         )
  )
))

# Convert the group column into a factor with the desired order
plot_data$Group <- factor(plot_data$Group, levels = c("Approved Targets", "Updated Approved", "Clinical Targets", "Updated Clinical", "Not Targets"))


ggplot(plot_data, aes(x = Group, y = Prediction_Score_10000m)) +
  geom_boxplot() +
  theme_minimal() +
  labs(
    x = NULL,
    y = "Prediction Scores"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  )


temp <- feature_data
temp$updated_only <- ifelse(temp$Label == 0 & temp$updated_approved == 1, 1, 0)

roc_curve <- roc(temp$updated_only, temp$Prediction_Score_100m)
auc(roc_curve)
plot(roc_curve)


temp <- plot_data[plot_data$Group == "Updated Approved", ]
