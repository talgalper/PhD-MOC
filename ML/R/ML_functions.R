### functions associated with ML work ###

# formats TTD data into positive set and negative pool. Only works from this directory on the TTD_data.RData object.
data_sets_from_TTD <- function(feature_matrix, ensembl, filter_terms = c("cancer", "carcinoma", "leukemia", "leukaemia", "neoplasm", "metastases", "tumour", "adenoma")) {
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
  all_target_genes <- strsplit(all_target_genes, "/")
  all_target_genes <- unlist(all_target_genes)
  all_target_genes <- unique(all_target_genes)
  all_target_genes <- trimws(all_target_genes)
  all_target_genes <- all_target_genes[all_target_genes != ""]

  # convert gene symbols into uniprot IDs
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
  pattern <- paste(filter_terms, collapse = "|")
  TTD_approved_cancer <- TTD_approved[grepl(pattern, TTD_approved$INDICATION, ignore.case = TRUE), ]
  TTD_approved_cancer <- TTD_approved_cancer[!is.na(TTD_approved_cancer$INDICATION), ]
  TTD_approved_cancer <- TTD_approved_cancer[grepl(c("small molecul"), TTD_approved_cancer$DRUGTYPE, ignore.case = TRUE), ] # molecule spelt wrong on purpose
  TTD_approved_cancer <- TTD_approved_cancer[!is.na(TTD_approved_cancer$DRUGTYPE), ]
  
  TTD_approved_cancer <- unique(TTD_approved_cancer$uniprotswissprot)
  TTD_approved_cancer <- TTD_approved_cancer[TTD_approved_cancer != "" & !is.na(TTD_approved_cancer)]
  
  feature_matrix$approved <- as.factor(ifelse(feature_matrix$Protein %in% TTD_approved_cancer, 1, 0))
  
  
  TTD_clinical <- merge(unique(TTD_master[, c(1,2,4,10,11)]), TTD_clinical, by = "TARGETID", all.y = T)
  pattern <- paste(filter_terms, collapse = "|")
  TTD_clinical_cancer <- TTD_clinical[grepl(pattern, TTD_clinical$INDICATION, ignore.case = TRUE), ]
  TTD_clinical_cancer <- TTD_clinical_cancer[!is.na(TTD_clinical_cancer$INDICATION), ]
  #looks like they did not use small molecule clinical trial drugs for validation, im doing it though
  TTD_clinical_cancer <- TTD_clinical_cancer[grepl(c("small molecul"), TTD_clinical_cancer$DRUGTYPE, ignore.case = TRUE), ] # molecule spelt wrong on purpose
  TTD_clinical_cancer <- TTD_clinical_cancer[!is.na(TTD_clinical_cancer$DRUGTYPE), ]
  
  TTD_clinical_cancer <- unique(TTD_clinical_cancer$uniprotswissprot)
  TTD_clinical_cancer <- TTD_clinical_cancer[TTD_clinical_cancer != "" & !is.na(TTD_clinical_cancer)]
  
  feature_matrix$clinical <- as.factor(ifelse(feature_matrix$Protein %in% TTD_clinical_cancer, 1, 0))
  
  # Split data into positive (approved drug targets) and negative (non-drug targets across any indication)
  positive_set <- feature_matrix[feature_matrix$approved == 1, ]
  negative_pool <- feature_matrix[feature_matrix$approved == 0 & !feature_matrix$Protein %in% TTD_approved$uniprotswissprot, ]
  
  return(list(positive_set = positive_set, negative_pool = negative_pool, feature_matrix = feature_matrix))
}


# Perform Random Forest training and test including parameters for tuning and parallelisation
RF_bagging <- function(feature_matrix, positive_set, negative_pool, 
                       ntrees = 1000, n_models = 100, 
                       track_iterations = TRUE, 
                       model_data_output_dir = NULL, 
                       parallel = FALSE, 
                       tuning = FALSE, 
                       set_mtry_tuning_grid = seq(2, 10, by = 1)) {
  
  suppressMessages({
    library(progressr)
    library(progress)
    library(foreach)
    library(doParallel)
    library(caret)
    library(randomForest)
    library(pROC)
    library(cli)
    library(doSNOW)
  })
  
  # Let the user know what is happening
  if (!parallel & tuning) {
    cat("WARNING: Tuning mtry without parallel threads will take a very LONG time...\n")
  } else if (parallel & !tuning) {
    cat("Running", n_models, "models using", detectCores() - 1, "threads, no tuning...\n")
  } else if (parallel & tuning) {
    cat("Running", n_models, "models with tuning, using", detectCores() - 1, "parallel threads...\n")
  } else {
    cat("Running", n_models, "models without tuning or parallel threads...\n")
  }
  
  # Prepare storage for predictions and importance
  predictions <- matrix(0, nrow = nrow(feature_matrix), ncol = n_models)
  n_features <- ncol(feature_matrix) - 3  
  importance_scores <- matrix(0, nrow = n_features, ncol = n_models)
  mtrys <- c()  # to store selected mtry values if tuning is used
  
  ##############################################################################
  #                              PARALLEL BRANCH                               #
  ##############################################################################
  if (isTRUE(parallel)) {
    # Create and register a cluster
    cl <- makeCluster(detectCores() - 1)
    registerDoSNOW(cl)
    
    # set up progress bar and allow to be used in foreach using doSNOW
    pb <- progress_bar$new(
      format = "[:bar] :current/:total (:percent) eta: :eta", 
      total  = n_models)
    progress <- function(n){
      pb$tick(tokens = list(model = rep(1:n_models)))
    }
    opts <- list(progress = progress)
    
    # Parallel loop with foreach
    # .combine='list' means each iteration returns a list, and we get a list-of-lists
    results_list <- foreach(
      i = seq_len(n_models),
      .combine = 'list',
      .packages = c("randomForest", "caret", "pROC"),
      .options.snow = opts
    ) %dopar% {
      set.seed(i)  # A different seed each iteration (partially reproducible)
      
      # Resample negatives for each iteration
      negative_set <- negative_pool[sample(seq_len(nrow(negative_pool)), 
                                           nrow(positive_set), replace = FALSE), ]
      training_data <- rbind(positive_set, negative_set)
      training_labels <- training_data$approved
      training_features <- training_data[, !names(training_data) %in% c("approved", "Protein", "clinical")]
      
      if (isTRUE(tuning)) {
        # Allow parallel CV
        tune_control <- trainControl(
          method = "cv", 
          number = 5, 
          verboseIter = FALSE, 
          allowParallel = FALSE
        )
        tune_grid <- expand.grid(mtry = set_mtry_tuning_grid)
        
        rf_tuned <- caret::train(
          x = training_features,
          y = as.factor(training_labels),
          method = "rf",
          metric = "Accuracy",
          tuneGrid = tune_grid,
          trControl = tune_control,
          ntree = ntrees
        )
        
        best_mtry <- rf_tuned$bestTune$mtry
        mtry_used <- best_mtry
        rf_model <- randomForest(
          x = training_features,
          y = as.factor(training_labels),
          mtry = best_mtry,
          ntree = ntrees,
          importance = TRUE
        )
      } else {
        # Standard randomForest
        mtry_used <- round(sqrt(ncol(training_features)))
        rf_model <- randomForest(
          x = training_features,
          y = as.factor(training_labels),
          mtry = mtry_used,
          ntree = ntrees,
          importance = TRUE
        )
      }
      
      # Predict on the entire feature_matrix
      pred_col <- predict(rf_model, 
                          feature_matrix[, !names(feature_matrix) %in% c("approved", "Protein", "clinical")], 
                          type = "prob")[, 2]
      
      # Extract feature importance
      imp_vals <- importance(rf_model)[, "MeanDecreaseGini"]
      
      # Return everything we need for combining later
      list(
        prediction = pred_col,
        importance = imp_vals,
        mtry_used = mtry_used
      )
    }
    close(pb)
    stopCluster(cl)
    
    # Collect results into predictions/importance_scores
    for (i in seq_len(n_models)) {
      predictions[, i] <- results_list[[i]]$prediction
      importance_scores[, i] <- results_list[[i]]$importance
      mtrys[i] <- results_list[[i]]$mtry_used
    }
    
    ##############################################################################
    #                              NON-PARALLEL BRANCH                           #
    ##############################################################################
  } else {
    # Simple progress bar (non-parallel)
    pb <- progress_bar$new(
      format = "[:bar] :current/:total (:percent) eta: :eta", 
      total  = n_models
    )
    
    for (i in seq_len(n_models)) {
      pb$tick()
      set.seed(i)
      
      # Resample negatives
      negative_set <- negative_pool[sample(seq_len(nrow(negative_pool)), 
                                           nrow(positive_set), replace = FALSE), ]
      training_data <- rbind(positive_set, negative_set)
      training_labels <- training_data$approved
      training_features <- training_data[, !names(training_data) %in% c("approved", "Protein", "clinical")]
      
      if (isTRUE(tuning)) {
        # Tuning
        tune_control <- trainControl(method = "cv", number = 5, verboseIter = FALSE, allowParallel = FALSE)
        tune_grid    <- expand.grid(mtry = set_mtry_tuning_grid)
        
        rf_tuned <- train(
          x = training_features,
          y = as.factor(training_labels),
          method = "rf",
          metric = "Accuracy",
          tuneGrid = tune_grid,
          trControl = tune_control,
          ntree = ntrees
        )
        
        best_mtry <- rf_tuned$bestTune$mtry
        mtrys[i]  <- best_mtry
        
        rf_model <- randomForest(
          x = training_features,
          y = as.factor(training_labels),
          mtry = best_mtry,
          ntree = ntrees,
          importance = TRUE
        )
      } else {
        # No tuning
        mtry_val   <- round(sqrt(ncol(training_features)))
        mtrys[i]   <- mtry_val
        
        rf_model <- randomForest(
          x = training_features,
          y = as.factor(training_labels),
          mtry = mtry_val,
          ntree = ntrees,
          importance = TRUE
        )
      }
      
      predictions[, i] <- predict(rf_model, 
                                  feature_matrix[, !names(feature_matrix) %in% c("approved", "Protein", "clinical")], 
                                  type = "prob")[, 2]
      imp_vals <- randomForest::importance(rf_model)[, "MeanDecreaseGini"]
      importance_scores[, i] <- imp_vals
    }
  }
  
  ##############################################################################
  #                 Summarise into feature_matrix + track iterations           #
  ##############################################################################
  if (isTRUE(track_iterations)) {
    # Which iterations do we track?
    track_iterations_list <- unique(c(1, 10^(0:floor(log10(n_models))), n_models))
    
    # Feature importance data frame
    feature_names <- colnames(feature_matrix)[
      !colnames(feature_matrix) %in% c("Protein", "approved", "clinical")
    ]
    importance_df <- data.frame(Features = feature_names)
    
    for (iter in track_iterations_list) {
      if (iter == 1) {
        importance_df[[paste0("Importance_", iter, "m")]] <- importance_scores[, 1]
      } else {
        importance_df[[paste0("Importance_", iter, "m")]] <- rowMeans(importance_scores[, seq_len(iter)])
      }
    }
    
    # Add columns for tracked predictions
    for (iter in track_iterations_list) {
      if (iter == 1) {
        feature_matrix[[paste0("Prediction_Score_", iter, "m")]] <- predictions[, 1]
      } else {
        feature_matrix[[paste0("Prediction_Score_", iter, "m")]] <- 
          rowMeans(predictions[, seq_len(iter)])
      }
    }
    
    # Evaluate performance predicting 'clinical' via ROC/AUC
    AUC <- data.frame(N_iterations = track_iterations_list)
    for (i in seq_along(track_iterations_list)) {
      iter <- track_iterations_list[i]
      roc_curve <- suppressMessages({
        roc(feature_matrix$clinical, feature_matrix[[paste0("Prediction_Score_", iter, "m")]])
      })
      AUC[i, "AUC"] <- pROC::auc(roc_curve)
    }
    
  } else {
    # Single final columns
    feature_matrix[[paste0("Prediction_Score_", n_models, "m")]] <- rowMeans(predictions)
    
    feature_names <- colnames(feature_matrix)[
      !colnames(feature_matrix) %in% c("Protein", "approved", "clinical")
    ]
    importance_df <- data.frame(Features = feature_names)
    importance_df[[paste0("Importance_", n_models, "m")]] <- rowMeans(importance_scores)
    
    # Evaluate performance
    roc_curve <- suppressMessages({
      roc(feature_matrix$clinical, feature_matrix[[paste0("Prediction_Score_", n_models, "m")]])
    })
    AUC <- pROC::auc(roc_curve)
  }
  
  ##############################################################################
  #                           (Optional) Save to .RData                         #
  ##############################################################################
  if (!is.null(model_data_output_dir)) {
    save(predictions, importance_scores, importance_df, AUC, feature_matrix, 
         file = model_data_output_dir)
  }
  
  ##############################################################################
  #                               Return results                                #
  ##############################################################################
  return(list(
    predictions = predictions,
    importance_scores = importance_scores,
    importance_df = importance_df,
    AUC = AUC,
    mtrys = mtrys
  ))
}

