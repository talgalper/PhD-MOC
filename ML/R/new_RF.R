library(randomForest)
library(caret)
library(progress)
library(pROC)

druggability <- read.csv("../Druggability_analysis/data_general/druggability_scores_annot.csv")
load("RData/full_fpocket_results.RData")
STRING_hsa <- read.csv("../Hierarchical_HotNet/STRING_data/STRING_physical_geneSymbol.csv")

results_master <- results_master[,-23]
results_master$Score <- (results_master$Score - min(results_master$Score)) / (max(results_master$Score) - min(results_master$Score))
results_master <- results_master[order(-results_master$`Druggability Score`, -results_master$Score), ]
results_master <- results_master[!duplicated(results_master$uniprot_id), ]

druggability <- subset(druggability, select = c("uniprot_gn_id", "CP_score", "highest_score"))
druggability <- unique(druggability)

feature_matrix <- merge(druggability, results_master[, -c(1,22)], by.x = "uniprot_gn_id", by.y = "uniprot_id", all = T)


library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

filters <- c("uniprotswissprot", "uniprot_gn_id", "uniprotsptrembl")
temp <- data.frame(uniprot_gn_id = feature_matrix$uniprot_gn_id)
for (filter in filters) {
  cat("Processing in BiomaRt using", filter, "as filter", "\n")
  
  gene_ids <- getBM(
    attributes = c(filter, "hgnc_symbol", "external_gene_name", "uniprot_gn_symbol"),
    filters = filter,
    values = feature_matrix$uniprot_gn_id,
    mart = ensembl)
  
  # Create the combined column in base R
  gene_ids$combined <- ifelse(
    gene_ids$hgnc_symbol != "" & !is.na(gene_ids$hgnc_symbol), 
    gene_ids$hgnc_symbol,
    ifelse(
      gene_ids$external_gene_name != "" & !is.na(gene_ids$external_gene_name),
      gene_ids$external_gene_name,
      ifelse(
        gene_ids$uniprot_gn_symbol != "" & !is.na(gene_ids$uniprot_gn_symbol),
        gene_ids$uniprot_gn_symbol,
        NA_character_
      )
    )
  )
  
  gene_ids <- gene_ids[, c(1,5)]
  colnames(gene_ids)[2] <- paste0("gene_symbol_", filter)
  gene_ids <- merge(feature_matrix, gene_ids, by.x = "uniprot_gn_id", by.y = filter, all.x = T)
  gene_ids <- gene_ids[, c(1,24)]
  gene_ids <- unique(gene_ids)
  
  temp <- merge(temp, gene_ids, by = "uniprot_gn_id")
}
rm(gene_ids, filter)


# Create the combined column in base R
temp$combined <- ifelse(
  temp$gene_symbol_uniprotswissprot != "" & !is.na(temp$gene_symbol_uniprotswissprot), 
  temp$gene_symbol_uniprotswissprot,
  ifelse(
    temp$gene_symbol_uniprot_gn_id != "" & !is.na(temp$gene_symbol_uniprot_gn_id),
    temp$gene_symbol_uniprot_gn_id,
    ifelse(
      temp$gene_symbol_uniprotsptrembl != "" & !is.na(temp$gene_symbol_uniprotsptrembl),
      temp$gene_symbol_uniprotsptrembl,
      NA_character_
    )
  )
)

gene_ids <- temp[, c(1,5)]
colnames(gene_ids)[2] <- "gene_symbol"
gene_ids <- merge(feature_matrix, gene_ids, by = "uniprot_gn_id", all.x = T)
gene_ids <- gene_ids[, c(1,24)]
gene_ids <- unique(gene_ids)

temp2 <- gene_ids[is.na(gene_ids$gene_symbol), ]

feature_matrix <- merge(gene_ids, feature_matrix, by = "uniprot_gn_id", all.y = T)

library(igraph)
STRING_net <- graph_from_data_frame(STRING_hsa, directed = F)

string_df <- data.frame(ENSG = V(STRING_net)$name,
                        degree = degree(STRING_net),
                        betweenness = betweenness(STRING_net),
                        closeness = closeness(STRING_net),
                        eigen_centrality = eigen_centrality(STRING_net)$vector,
                        page_rank = page_rank(STRING_net)$vector)

feature_matrix <- merge(feature_matrix, string_df, by.x = "gene_symbol", by.y = "ENSG", all = T)

feature_matrix <- feature_matrix[complete.cases(feature_matrix), ]
feature_matrix <- unique(feature_matrix)

feature_matrix[c(24:28)] <- lapply(feature_matrix[c(24:28)], function(x) {
  (x - min(x)) / (max(x) - min(x))
})

feature_matrix <- feature_matrix[order(-feature_matrix$`Druggability Score`, -feature_matrix$Score), ]



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


TTD_approved <- merge(unique(TTD_master[, c(1,2,3,9)]), TTD_approved, by = "TARGETID", all.y = T)
cancer_terms <- c("cancer", "carcinoma", "leukemia", "neoplasm", "metastases", "tumour")
pattern <- paste(cancer_terms, collapse = "|")
TTD_approved_cancer <- TTD_approved[grepl(pattern, TTD_approved$INDICATION, ignore.case = TRUE), ]
TTD_approved_cancer <- TTD_approved_cancer[!is.na(TTD_approved_cancer$INDICATION), ]

TTD_approved_cancer <- unique(TTD_approved_cancer$GENENAME)
TTD_approved_cancer <- TTD_approved_cancer[TTD_approved_cancer != "" & !is.na(TTD_approved_cancer)]

feature_matrix$approved <- as.factor(ifelse(feature_matrix$gene_symbol %in% TTD_approved_cancer, 1, 0))


TTD_clinical <- merge(unique(TTD_master[, c(1,2,3,9)]), TTD_clinical, by = "TARGETID", all.y = T)
cancer_terms <- c("cancer", "carcinoma", "leukemia", "neoplasm", "metastases")
pattern <- paste(cancer_terms, collapse = "|")
TTD_clinical_cancer <- TTD_clinical[grepl(pattern, TTD_clinical$INDICATION, ignore.case = TRUE), ]
TTD_clinical_cancer <- TTD_clinical_cancer[!is.na(TTD_clinical_cancer$INDICATION), ]

TTD_clinical_cancer <- unique(TTD_clinical_cancer$GENENAME)
TTD_clinical_cancer <- TTD_clinical_cancer[TTD_clinical_cancer != "" & !is.na(TTD_clinical_cancer)]

feature_matrix$clinical <- as.factor(ifelse(feature_matrix$gene_symbol %in% TTD_clinical_cancer, 1, 0))


# Split data into positive (approved drug targets) and negative (non-drug targets across any indication)
positive_set <- feature_matrix[feature_matrix$approved == 1, ]
negative_pool <- feature_matrix[!feature_matrix$gene_symbol %in% all_target_genes, ]





library(corrplot)

# Calculate correlations between features
# Exclude non-numeric columns like 'Protein', 'Label', or 'validation'
numeric_features <- feature_matrix[, sapply(feature_matrix, is.numeric)]
correlation_matrix <- cor(numeric_features, use = "pairwise.complete.obs")

# Plot the correlation matrix
corrplot(correlation_matrix, method = "color",
         tl.cex = 0.7, tl.col = "black", tl.srt = 45, 
         mar = c(0, 0, 0, 0))





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
  training_features <- training_data[, !names(training_data) %in% c("gene_symbol", "uniprot_gn_id", "approved", "clinical")]
  
  # Train Random Forest
  rf_model <- randomForest(
    x = training_features,
    y = as.factor(training_labels),
    mtry = round(sqrt(ncol(training_features))),
    ntree = ntrees,
    importance = TRUE)
  
  # Predict on the entire dataset
  predictions[, i] <- predict(rf_model, feature_matrix[, !names(feature_matrix) %in% c("approved", "Protein", "clinical")], type = "prob")[, 2]
  
  pb$tick()
}


# Average predictions across all models
final_predictions <- rowMeans(predictions)
# Attach predictions to the dataset
feature_matrix$Prediction_Score <- final_predictions

model_scores <- subset(feature_matrix, select = c("gene_symbol", "Prediction_Score", "approved", "clinical"))
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


save(roc_curve, importance_df, model_scores, feature_matrix, positive_set, negative_pool, file = "RData/new_RF_model.RData")
load("RData/new_RF_model.RData")


plot_data <- subset(feature_matrix, select = c("gene_symbol", "Prediction_Score", "approved", "clinical"))


plot_data$Group <- with(plot_data, ifelse(
  approved == 1, "Approved Targets",
  ifelse(clinical == 1, "Clinical Targets", "Not Targets")
))


# Convert the group column into a factor with the desired order
plot_data$Group <- factor(plot_data$Group, levels = c("Approved Targets", "Clinical Targets", "Not Targets"))


ggplot(plot_data, aes(x = Group, y = Prediction_Score)) +
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
