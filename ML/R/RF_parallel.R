RF_bagging <- function(feature_matrix, positive_set, negative_pool, ntrees = 1000, n_models = 100, track_iterations = TRUE, model_data_output_dir = NULL) {
  
  suppressMessages({
    library(progressr)
    library(foreach)
    library(doParallel)
    library(caret)
    library(randomForest)
    library(pROC)
    library(cli)
  })
  
  handlers("cli") # Set up progress handlers
  cl <- makeCluster(detectCores() - 1)
  registerDoParallel(cl)
  
  predictions <- matrix(0, nrow = nrow(feature_matrix), ncol = n_models)  # For storing predictions
  n_features <- ncol(feature_matrix) - 3  # For storing importance scores
  importance_scores <- matrix(0, nrow = n_features, ncol = n_models)
  
  tune_grid <- expand.grid(
    mtry = seq(2, 10, by = 1)  # Only include 'mtry'
  )
  
  mtrys <- c()
  
  # Use progressr to track the loop progress
  with_progress({
    p <- progressor(along = 1:n_models) # Create a progressor
    
    # Train multiple RF models using random negative samples
    for (i in 1:n_models) {
      p(sprintf("Running model %d/%d", i, n_models)) # Update progress
      
      set.seed(i)  # Different seed for each iteration
      
      # Resample negatives for each iteration
      negative_set <- negative_pool[sample(1:nrow(negative_pool), nrow(positive_set), replace = FALSE), ]
      training_data <- rbind(positive_set, negative_set)
      training_labels <- training_data$approved
      training_features <- training_data[, !names(training_data) %in% c("approved", "Protein", "clinical")]
      
      # Perform hyperparameter tuning
      tune_control <- trainControl(method = "cv", number = 5, verboseIter = FALSE)
      
      rf_tuned <- train(
        x = training_features,
        y = as.factor(training_labels),
        method = "rf",
        metric = "Accuracy",
        tuneGrid = tune_grid,  # Only 'mtry'
        trControl = tune_control,
        ntree = ntrees  # Pass 'ntree' directly here
      )
      
      # Use best tuned parameters
      best_mtry <- rf_tuned$bestTune$mtry
      mtrys <- append(mtrys, best_mtry)

      # Train Random Forest with tuned parameters
      rf_model <- randomForest(
        x = training_features,
        y = as.factor(training_labels),
        mtry = best_mtry,
        ntree = ntrees,
        importance = TRUE
      )
      
      # Predict on the entire dataset
      predictions[, i] <- predict(rf_model, feature_matrix[, !names(feature_matrix) %in% c("approved", "Protein", "clinical")], type = "prob")[, 2]
      
      # Extract feature importance
      importance <- importance(rf_model)
      importance_scores[, i] <- importance[, "MeanDecreaseGini"]
    }
  })
  
  # Remaining part of the function unchanged...
  
  # Attach tracked predictions to the dataset
  if (isTRUE(track_iterations)) {
    track_iterations <- unique(c(1, 10^(0:floor(log10(n_models))), n_models))
    
    # Calculate mean importance across rows for tracked iterations
    importance_df <- data.frame(Features = colnames(feature_matrix)[!colnames(feature_matrix) %in% c("Protein", "approved", "clinical")])
    for (iter in track_iterations) {
      if (iter == 1) {
        importance_df[[paste0("Importance_", iter, "m")]] <- importance_scores[, 1]
      } else {
        importance_df[[paste0("Importance_", iter, "m")]] <- rowMeans(importance_scores[, 1:iter])
      }
    }
    
    # Calculate mean prediction score across rows for tracked iterations
    for (iter in track_iterations) {
      if (iter == 1) {
        feature_matrix[[paste0("Prediction_Score_", iter, "m")]] <- predictions[, 1]
      } else {
        feature_matrix[[paste0("Prediction_Score_", iter, "m")]] <- rowMeans(predictions[, 1:iter])
      }
    }
    
    # Evaluate performance predicting clinical oncology small molecule drugs
    AUC <- data.frame(N_iterations = c(track_iterations))
    for (i in seq_along(track_iterations)) {
      iter <- track_iterations[i]
      roc_curve <- suppressMessages({roc(feature_matrix$clinical, feature_matrix[[paste0("Prediction_Score_", iter, "m")]])})
      AUC[i, "AUC"] <- auc(roc_curve)
    }
  } else {
    feature_matrix[[paste0("Prediction_Score_", n_models, "m")]] <- rowMeans(predictions)
    importance_df <- data.frame(Features = colnames(feature_matrix)[!colnames(feature_matrix) %in% c("Protein", "approved", "clinical")])
    importance_df[[paste0("Importance_", n_models, "m")]] <- rowMeans(importance_scores)
    roc_curve <- roc(feature_matrix$clinical, feature_matrix[[paste0("Prediction_Score_", n_models, "m")]])
    AUC <- auc(roc_curve)
  }
  
  # Save RData object (optional)
  if (!is.null(model_data_output_dir)) {
    save(predictions, importance_scores, importance_df, AUC, feature_matrix, file = model_data_output_dir)
  }
  
  return(list(predictions = predictions, importance_scores = importance_scores, importance_df = importance_df, AUC = AUC, mtrys = mtrys))
  gc()
  stopCluster(cl)
}

library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
feature_matrix <- read.table("data/feature_matrix.txt", sep = "\t", header = T)
training_data <- data_sets_from_TTD(ensembl)

RF_results2 <- RF_bagging(training_data$feature_matrix, 
                         positive_set = training_data$positive_set, negative_pool = training_data$negative_pool, 
                         n_models = 1000, 
                         track_iterations = T)

