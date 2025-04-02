##############################################################################
#            run ML bagging function on dezso matrix only                    #
##############################################################################
library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
feature_matrix <- read.table("data/feature_matrix.txt", sep = "\t", header = T)
training_data <- data_sets_from_TTD(feature_matrix, ensembl)

start <- Sys.time()
ML_bagging_results <- ML_bagging(feature_matrix = training_data$feature_matrix,
                                 positive_set = training_data$positive_set, 
                                 negative_pool = training_data$negative_pool,
                                 models = c("glmnet", "rf", "svmRadial", "knn", "nb", "nnet", "xgbTree"),
                                 n_models = 100,
                                 big_grid = FALSE)
print(Sys.time() - start)
rm(start)

save(ML_bagging_results, file = "~/OneDrive - RMIT University/PhD/large_git_files/ML/ML_bagging_100itr.RData")
load("~/OneDrive - RMIT University/PhD/large_git_files/ML/ML_bagging_100itr.RData")

# plot cross fold validation results
plot_result <- plot_results(ML_bagging_results = ML_bagging_results,
                            return_plot_data = TRUE)

# extract feature importance
feature_importance <- extract_feature_importance(ML_bagging_results, plot_feature_importance = TRUE)

# get model performance metrics on predicting clinical targets
model_prediction_results <- model_metrics(feature_matrix = feature_matrix, 
                                          ML_bagging_results = ML_bagging_results,
                                          plot_results = TRUE)



###################################################################################
# run ML bagging function on dezso matrix with druggability and new STRING scores #
#               reduced number of models used to decrease run time                 #
###################################################################################
# append new data to feature matrix
library(readxl)
feature_data <- read.table("data/feature_matrix.txt", sep = "\t", header = T)
approved_targets <- read_xls("data/approved_targets.xls", sheet = 1)
clinical_targets <- read_xls("data/approved_targets.xls", sheet = 2)

druggability <- read.csv("../Druggability_analysis/data_general/druggability_scores_annot.csv")
load("RData/full_fpocket_results.RData")
load("RData/human_string_PPI_metrics.RData")
load("RData/human_STRING_top10pct_PPI_metrics.RData")

results_master <- results_master[,-23]
results_master$Score <- (results_master$Score - min(results_master$Score)) / (max(results_master$Score) - min(results_master$Score))
results_master <- results_master[order(-results_master$`Druggability Score`, -results_master$Score), ]
results_master <- results_master[!duplicated(results_master$uniprot_id), ]
colnames(results_master)[3:21] <- paste0("pocket_", colnames(results_master)[3:21])

druggability <- subset(druggability, select = c("uniprot_gn_id", "CP_score", "highest_score"))
druggability <- unique(druggability)

colnames(STRING_top10pct_df) <- paste0(colnames(STRING_top10pct_df), "_10pct")

myFeatures <- merge(druggability, results_master[, -c(1,22)], by.x = "uniprot_gn_id", by.y = "uniprot_id", all = T)
myFeatures <- merge(myFeatures, STRING_top10pct_df, by.x = "uniprot_gn_id", by.y = "ENSG_10pct", all = T)
myFeatures <- merge(myFeatures, string_df, by.x = "uniprot_gn_id", by.y = "ENSG", all = T)
myFeatures[is.na(myFeatures)] <- 0

myFeatures[23:ncol(myFeatures)] <- lapply(myFeatures[23:ncol(myFeatures)], function(x) {
  (x - min(x)) / (max(x) - min(x))
})

# correlation matrix between similar features and updated data
library(corrplot)
library(reshape2)
correlation_matrix <- merge(feature_data, myFeatures, by.x = "Protein", by.y = "uniprot_gn_id", all.x = T)
correlation_matrix <- correlation_matrix[, sapply(correlation_matrix, is.numeric)]
correlation_matrix <- cor(correlation_matrix, use = "pairwise.complete.obs")
correlation_matrix <- correlation_matrix[-c(1:65), ]

corrplot(correlation_matrix, method = "color",
         tl.cex = 0.7, tl.col = "black", tl.srt = 45)

corr_melt <- melt(correlation_matrix)
corr_melt <- corr_melt[order(-corr_melt$value), ]
colnames(corr_melt) <- c("myFeatures", "features", "pearson_corr")

feature_data$approved <- as.factor(ifelse(feature_data$Protein %in% approved_targets$Protein, 1, 0))
feature_data$clinical <- as.factor(ifelse(feature_data$Protein %in% clinical_targets$Protein, 1, 0))

updated_feature_data <- merge(feature_data, myFeatures, by.x = "Protein", by.y = "uniprot_gn_id", all.x = T)
updated_feature_data[is.na(updated_feature_data)] <- 0

positive_set <- updated_feature_data[updated_feature_data$approved == 1, ]

library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
training_data <- data_sets_from_TTD(updated_feature_data, ensembl)

# run bagging
# use original positive training set and clinical with updated negative pool
start <- Sys.time()
ML_bagging_results <- ML_bagging(feature_matrix = updated_feature_data,
                                 positive_set = positive_set, 
                                 negative_pool = training_data$negative_pool, 
                                 n_models = 1000,
                                 models = c("glmnet", "rf", "svmRadial", "nnet"))
print(Sys.time() - start)
rm(start)


# plot cross fold validation results
plot_result <- plot_results(ML_bagging_results = ML_bagging_results,
                            return_plot_data = TRUE)

# extract feature importance
feature_importance <- extract_feature_importance(ML_bagging_results, plot_feature_importance = TRUE)

# get model performance metrics on predicting clinical targets
model_prediction_results <- model_metrics(feature_matrix = updated_feature_data, 
                                          ML_bagging_results = ML_bagging_results,
                                          plot_results = TRUE)

save(ML_bagging_results, feature_importance, model_prediction_results, plot_result, updated_feature_data,
     file = "~/OneDrive - RMIT University/PhD/large_git_files/ML/ML_bagging_more_features(100itr).RData")
load("~/OneDrive - RMIT University/PhD/large_git_files/ML/ML_bagging_more_features(100itr).RData")



## compare results with network methods
#HHnet_RS <- read.csv("../BRCA_pipe/latest_run/OG_rank/HHnet_rank_sensitivity_top10.csv", row.names = 1)
HHnetEnrich_RS <- read.csv("../BRCA_pipe/latest_run/OG_rank/avg_RS_HHnetEnrich.csv")
PCSF_RS <- read.csv("../BRCA_pipe/latest_run/OG_rank/avg_RS_PCSF.csv")

feature_data_scores_appended <- model_prediction_results$feature_data_scores_appended


targets <- read.csv("../Druggability_analysis/data_general/target_all_dbs.csv")
targets <- targets[, c(2:4)]
targets <- targets[!duplicated(targets$ensembl_gene_id), ]

temp <- feature_data_scores_appended[feature_data_scores_appended$Protein %in% targets$Uniprot_ID, ]
temp <- temp[, c(1,105:108)]
temp <- merge(targets, temp, by.x = "Uniprot_ID", by.y = "Protein")

temp2 <- feature_data_scores_appended[, c(1,105:108)]
temp2 <- temp2[order(-temp2$Prediction_Score_rf), ]
temp2 <- temp2[temp2$Protein %in% targets$Uniprot_ID, ]
rank <- rownames(temp2)
temp2 <- merge(targets, temp2, by.x = "Uniprot_ID", by.y = "Protein")
temp2 <- temp2[!duplicated(temp2$drugBank_target), ]
temp2 <- temp2[order(-temp2$Prediction_Score_rf), ]
rownames(temp2) <- rank

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# uniprot_ids <- getBM(
#   attributes = c("hgnc_symbol", "uniprotswissprot", "uniprot_gn_id"),
#   filters = "hgnc_symbol",
#   values = HHnet_RS$external_gene_name,
#   mart = ensembl)
# 
# HHnet_RS <- merge(uniprot_ids, HHnet_RS, by.x = "hgnc_symbol", by.y = "external_gene_name", all.y = T)
# HHnet_RS <- merge(HHnet_RS, feature_data_scores_appended[, c(1,105:108)], by.x = "uniprot_gn_id", by.y = "Protein", all.y = )
# HHnet_RS <- unique(HHnet_RS)
# HHnet_RS <- HHnet_RS[order(-HHnet_RS$count), ]
# HHnet_RS <- HHnet_RS[!duplicated(HHnet_RS$hgnc_symbol) & HHnet_RS$uniprot_gn_id != "", ]


uniprot_ids <- getBM(
  attributes = c("external_gene_name", "uniprotswissprot", "uniprot_gn_id"),
  filters = "external_gene_name",
  values = HHnetEnrich_RS$gene,
  mart = ensembl)

HHnetEnrich_RS <- merge(uniprot_ids, HHnetEnrich_RS, by.x = "external_gene_name", by.y = "gene", all.y = T)
HHnetEnrich_RS <- merge(HHnetEnrich_RS, feature_data_scores_appended[, c(1,105:108)], by.x = "uniprot_gn_id", by.y = "Protein", all.y = )
HHnetEnrich_RS <- unique(HHnetEnrich_RS)
HHnetEnrich_RS <- HHnetEnrich_RS[order(HHnetEnrich_RS$avg_rank), ]
HHnetEnrich_RS <- HHnetEnrich_RS[!duplicated(HHnetEnrich_RS$external_gene_name) & HHnetEnrich_RS$uniprot_gn_id != "", ]

uniprot_ids <- getBM(
  attributes = c("external_gene_name", "uniprotswissprot", "uniprot_gn_id"),
  filters = "external_gene_name",
  values = PCSF_RS$gene,
  mart = ensembl)

PCSF_RS <- merge(uniprot_ids, PCSF_RS, by.x = "external_gene_name", by.y = "gene", all.y = T)
PCSF_RS <- merge(PCSF_RS, feature_data_scores_appended[, c(1,105:108)], by.x = "uniprot_gn_id", by.y = "Protein", all.y = )
PCSF_RS <- unique(PCSF_RS)
PCSF_RS <- PCSF_RS[order(PCSF_RS$avg_rank), ]
PCSF_RS <- PCSF_RS[!duplicated(PCSF_RS$external_gene_name) & PCSF_RS$uniprot_gn_id != "", ]

table(model_prediction_results$feature_data_scores_appended$Prediction_Score_rf >= 0.5)
table(model_prediction_results$feature_data_scores_appended$Prediction_Score_rf >= 0.5 & 
        model_prediction_results$feature_data_scores_appended$clinical == 1)
table(model_prediction_results$feature_data_scores_appended$Prediction_Score_rf >= 0.5 & 
        model_prediction_results$feature_data_scores_appended$clinical != 1 &
        model_prediction_results$feature_data_scores_appended$approved != 1)


updated_feature_data[updated_feature_data$Protein == "P11511", "approved"]


predicted_targets <- model_prediction_results$feature_data_scores_appended$Protein[model_prediction_results$feature_data_scores_appended$Prediction_Score_rf >= 0.5]

venn_data <- list(#HHnet_RS = HHnet_RS$uniprot_gn_id,
                  HHnetEnrich = HHnetEnrich_RS$uniprot_gn_id,
                  PCSF = PCSF_RS$uniprot_gn_id,
                  `ML predicted` = predicted_targets)

library(venn)
library(RColorBrewer)
venn(venn_data, 
     ellipse = T, 
     zcolor = brewer.pal(length(venn_data), name = "Dark2"),
     box = FALSE,
     ilabels = "counts",
     sncs = 2,
     ilcs = 2)



# intersection of RS genes and ML predicted targets
rownames(HHnetEnrich_RS) <- NULL
temp3 <- HHnetEnrich_RS[HHnetEnrich_RS$uniprot_gn_id %in% predicted_targets, ]
temp3 <- temp3[order(temp3$avg_rank), ]
#temp3 <- temp3[!temp3$uniprot_gn_id %in% training_data$positive_set$Protein, ]
not_ml_pred <- HHnetEnrich_RS[!HHnetEnrich_RS$uniprot_gn_id %in% predicted_targets, ]

temp3 <- PCSF_RS[PCSF_RS$uniprot_gn_id %in% predicted_targets, ]
#temp3 <- temp3[!temp3$uniprot_gn_id %in% training_data$positive_set$Protein, ]

# indications for temp3
load("~/Desktop/temp.RData")
# TTD_master <- TTD_master[grep("breast", TTD_master$INDICATION, ignore.case = T), ]
# TTD_master <- TTD_master[TTD_master$HighestClinicalStatus == "Approved", ]
# TTD_master <- unique(TTD_master$all_target_genes)
# TTD_master <- na.omit(TTD_master)
# TTD_master <- strsplit(TTD_master, ";")
# TTD_master <- unlist(TTD_master)
# TTD_master <- strsplit(TTD_master, "/")
# TTD_master <- unlist(TTD_master)
# TTD_master <- unique(TTD_master)
# TTD_master <- trimws(TTD_master)

temp4 <- TTD_master[TTD_master$all_target_genes %in% temp3$external_gene_name, ]
temp4 <- temp4[!duplicated(temp4[, c("all_target_genes", "INDICATION")]), ]
length(unique(temp4$all_target_genes))

# get cancer, non-cancer and unknown related targets
pattern <- paste(c("cancer", "carcinoma", "leukemia", "leukaemia", "neoplasm", "metastases", "tumour", "adenoma", "sarcoma"), collapse = "|")
cancer_related <- temp4[grepl(pattern, temp4$INDICATION, ignore.case = TRUE), ]
length(unique((cancer_related$all_target_genes)))
nonCancer_related <- temp4[!temp4$all_target_genes %in% unique(cancer_related$all_target_genes), ]
unknown_indication <- nonCancer_related[is.na(nonCancer_related$INDICATION) & !duplicated(nonCancer_related$all_target_genes), ]
nonCancer_related <- nonCancer_related[!nonCancer_related$all_target_genes %in% unknown_indication$all_target_genes, ]
nonCancer_related <- temp3[temp3$external_gene_name %in% unique(nonCancer_related$all_target_genes), ]
unknown_indication <- temp3[temp3$external_gene_name %in% unique(unknown_indication$all_target_genes), ]
no_TTD <- temp3[!temp3$external_gene_name %in% unique(TTD_master$all_target_genes), ]
cancer_related <- temp3[temp3$external_gene_name %in% unique(cancer_related$all_target_genes), ]


# total breast cancer related drugs
breast_targets <- TTD_master[grepl("breast", TTD_master$INDICATION, ignore.case = TRUE), ]
breast_targets <- breast_targets[breast_targets$HighestClinicalStatus == "Approved", ]

breast_targets <- unique(breast_targets$all_target_genes)
breast_targets <- na.omit(breast_targets)
breast_targets <- strsplit(breast_targets, ";")
breast_targets <- unlist(breast_targets)
breast_targets <- strsplit(breast_targets, "/")
breast_targets <- unlist(breast_targets)
breast_targets <- unique(breast_targets)
breast_targets <- trimws(breast_targets)

table(breast_targets %in% cancer_related$hgnc_symbol)


rank <- rownames(not_ml_pred)
not_ml_pred <- merge.data.table(not_ml_pred, PubTator_counts, by.x = "hgnc_symbol", by.y = "symbol", all.x = T)
not_ml_pred <- not_ml_pred[order(not_ml_pred$avg_rank), ]
rownames(not_ml_pred) <- rank







