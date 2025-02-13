
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


ML_bagging <- function(n_models, feature_matrix, positive_set, negative_pool, models, big_grid = FALSE) {
  
  suppressMessages({
    library(caret)
    library(progress)
    library(foreach)
    library(doParallel)
    library(doSNOW)
    library(pROC)
  })
  
  # Define fixed numeric features based on feature_matrix
  numeric_features <- names(feature_matrix)[sapply(feature_matrix, is.numeric)]
  
  model_predictions <- list()
  for (model in models) {
    model_predictions[[model]] <- matrix(0, nrow = nrow(feature_matrix), ncol = n_models)
  }
  
  # Initialize resamples list
  resamples_list <- vector("list", n_models)
  
  # Set up parallel backend using doSNOW and 80% of cpu capacity
  cl <- makeSOCKcluster(min(round(detectCores()*0.8), n_models))
  registerDoSNOW(cl)
  on.exit(stopCluster(cl))
  
  
  # Export necessary variables and functions to the workers
  snow::clusterExport(cl,
                        list = c("numeric_features", "feature_matrix", "positive_set", "negative_pool", "models"),
                        envir = environment())
  
  # setup progress bar for parallelisation
  pb <- progress_bar$new(
    format = "[:bar] :current/:total (:percent) eta: :eta", 
    total  = n_models)
  progress <- function(n){
    pb$tick(tokens = list(model = rep(1:n_models)))
  }
  opts <- list(progress = progress)
  
  parallel_result <- foreach(
    iter = 1:n_models, 
    .options.snow = opts, 
    .combine = "c",
    .packages = c("caret", "pROC")
  ) %dopar% {
    
    set.seed(iter)
    
    negative_set <- negative_pool[sample(seq_len(nrow(negative_pool)), 
                                         nrow(positive_set), replace = FALSE), ]
    training_set <- rbind(positive_set, negative_set)
    training_labels <- training_set$approved
    training_labels <- factor(
      training_labels,   
      levels = c(0, 1),  
      labels = c("No", "Yes") 
    )
    training_features <- training_set[, !names(training_set) %in% c("approved", "Protein", "clinical")]
    training_features_numeric <- training_features[, numeric_features]
    
    
    # Train control with 10-fold cross-validation
    train_control <- trainControl(method = "cv", number = 10, 
                                  summaryFunction = multiClassSummary, 
                                  classProbs = TRUE, allowParallel = FALSE)
    
    n_features <- ncol(feature_matrix) - 3 # Adjust as per your data
    
    if (isTRUE(big_grid)) {
      tune_grids <- list(
        glmnet = expand.grid(
          alpha = c(0, 0.25, 0.5, 0.75, 1),
          lambda = 10^seq(-5, 2, length = 15)
        ),
        
        rf = expand.grid(
          mtry = seq(
            from = max(1, round(sqrt(n_features)) - 5), 
            to = min(n_features, round(sqrt(n_features)) + 5), 
            by = 1
          )
        ),
        
        svmRadial = expand.grid(
          sigma = c(0.001, 0.01, 0.1, 1, 10, 100),
          C = c(0.01, 0.1, 1, 10, 100, 1000)
        ),

        knn = expand.grid(
          k = seq(
            from = max(1, round(sqrt(n_features)) - 5),
            to   = round(sqrt(n_features)) + 5,
            by   = 2
          )
        ),
        
        nb = expand.grid(
          fL        = c(0, 0.25, 0.5, 1, 2),   # Laplace smoothing
          usekernel = c(TRUE, FALSE),
          adjust    = c(0.25, 0.5, 1, 2, 4)
        ),
        
        nnet = expand.grid(
          size  = c(5, 10, 15, 20, 30),
          decay = c(1e-4, 1e-3, 1e-2, 1e-1, 1)
        ),
        
        xgbTree = expand.grid(
          nrounds         = c(50, 100, 150, 200),
          eta             = c(0.001, 0.01, 0.1, 0.3, 0.5, 1),
          max_depth       = c(3, 5, 7, 9),
          gamma           = c(0, 1, 2, 5),
          colsample_bytree= c(0.5, 0.7, 0.8, 1),
          min_child_weight= c(1, 3, 5),
          subsample       = c(0.5, 0.75, 1)
        )
      )
    } else {
      # Hyperparameter grids
      tune_grids <- list(
        glmnet = expand.grid(alpha = c(0, 0.5, 1), lambda = 10^seq(-4, 1, length = 10)),
        rf = expand.grid(mtry = seq(round(sqrt(n_features)) - 3, round(sqrt(n_features)) + 3, by = 1)),
        svmRadial = expand.grid(sigma = c(0.01, 0.1, 1, 10), C = c(0.1, 1, 10, 100)),
        knn = expand.grid(k = seq(round(sqrt(n_features)) - 3, round(sqrt(n_features)) + 3, by = 1)),  
        nb = expand.grid(fL = c(0, 0.5, 1, 2), usekernel = c(TRUE, FALSE), adjust = c(0.5, 1, 2)),
        nnet = expand.grid(size = c(5, 10, 15), decay = c(0.001, 0.01, 0.1, 1)),
        xgbTree = expand.grid(
          nrounds = c(50, 100, 150),
          eta = c(0.01, 0.1, 1),
          max_depth = c(3, 6, 9),
          gamma = c(0, 1),
          colsample_bytree = c(0.5, 0.8, 1),
          min_child_weight = c(1, 5),
          subsample = c(1)
        )
      )
    }

    # initialise local data
    local_model_predictions <- list()
    # Train all models on this iteration's data within the silenced block
    this_iter_models <- list()
    for (model_name in models) {
      set.seed(iter)
      
      if (model_name == "xgbTree") {
        Sys.setenv(XGBOOST_LOG_LEVEL = "0")
        fit <- train(
          x = training_features_numeric, 
          y = training_labels, 
          method = model_name, 
          trControl = train_control, 
          tuneGrid = tune_grids[[model_name]],
          metric = "Accuracy", 
          verbose = FALSE,
          nthread = 1
        )
      } else {
        fit <- train(
          x = training_features_numeric, 
          y = training_labels, 
          method = model_name, 
          trControl = train_control, 
          tuneGrid = tune_grids[[model_name]],
          metric = "Accuracy"
        )
      }
      this_iter_models[[model_name]] <- fit
    }
    
    local_resamples <- resamples(this_iter_models)
    
    
    for (model_name in models) {
      # Predict probabilities
      predicted_probs <- predict(
        this_iter_models[[model_name]],
        newdata = feature_matrix[, numeric_features],
        type = "prob"
      )[, 2]
      
      # Assign predictions to the local model_predictions list
      local_model_predictions[[model_name]] <- predicted_probs
    }
    
    # Return predictions and resamples for this iteration to be combined
    results <- list(model_predictions = local_model_predictions,
                       resamples = local_resamples,
                       iter_models = this_iter_models)
    list(results)
  } # dopar end

  # Combine the foreach results into model_predictions and resamples
  iter_models <- list()
  for (iter in 1:n_models) {
    iteration_result <- parallel_result[[iter]]
    iter_models[[iter]] <- parallel_result[[iter]][["iter_models"]]
    
    # Assign predictions
    for (model_name in models) {
      model_predictions[[model_name]][, iter] <- iteration_result$model_predictions[[model_name]]
    }
    # Assign resamples
    resamples_list[[iter]] <- iteration_result$resamples
  }
  
  # mean prediction across iterations
  average_predictions <- lapply(model_predictions, function(pred_matrix) {
    rowMeans(pred_matrix)
  })
  
 return(list(iter_models = iter_models,
             model_predictions = model_predictions,
             average_predictions = average_predictions,
             resamples = resamples_list)
        )
} # function end


##############################################################################
#                       Plot cross validation results                        #
##############################################################################
plot_results <- function(ML_bagging_results, return_plot_data = FALSE) {
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(reshape2)
  
  # Initialize a list to store performance metrics
  performance_list <- list()
  
  # Loop through each iteration and extract all available metrics
  for (i in seq_along(ML_bagging_results$resamples)) {
    resample <- ML_bagging_results$resamples[[i]]
    all_metrics <- summary(resample)$statistics
    
    for (metric in names(all_metrics)) {
      metric_values <- all_metrics[[metric]][, "Mean"]
      
      temp_df <- data.frame(
        Iteration = i,
        Model = names(metric_values),
        Metric = metric,
        Value = as.numeric(metric_values),
        stringsAsFactors = FALSE
      )
      
      performance_list[[length(performance_list) + 1]] <- temp_df
    }
  }
  
  # Combine all performance data
  performance_df <- bind_rows(performance_list)
  
  # Compute summary statistics (Mean & SD) for each metric across models
  summary_stats <- performance_df %>%
    group_by(Model, Metric) %>%
    summarise(
      Mean = mean(Value, na.rm = TRUE),
      SD = sd(Value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    as.data.frame()
  
  ### ANALYZE HYPERPARAMETERS ###
  hyperparameters_list <- list()
  
  for (i in seq_along(ML_bagging_results$iter_models)) {
    iter_models <- ML_bagging_results$iter_models[[i]]
    
    for (model_name in names(iter_models)) {
      model <- iter_models[[model_name]]
      best_tune <- model$bestTune
      
      if (!is.null(best_tune)) {
        best_tune$Iteration <- i
        best_tune$Model <- model_name
        hyperparameters_list[[length(hyperparameters_list) + 1]] <- best_tune
      }
    }
  }
  
  # Base R Boxplot of Accuracy values
  print(
    ggplot(performance_df %>% filter(Metric == "Accuracy"), aes(x = Model, y = Value)) +
    geom_boxplot(fill = "steelblue") +
    labs(title = "Distribution of Accuracy Across Iterations", y = "Accuracy") +
    theme(axis.text = element_text(size = 20, colour = "black"),
          axis.title.x = element_text(size = 22, face = "bold", margin = margin(t=10)),
          axis.title.y = element_text(size = 22, face = "bold", margin = margin(r=10)),
          panel.grid = element_blank(),
          panel.background = element_blank()))
  
  # ggplot2 Line Plot
  print(
    ggplot(performance_df %>% filter(Metric == "Accuracy"), aes(x = Iteration, y = Value, color = Model)) +
    geom_line(alpha = 0.6) +
    geom_smooth(method = "loess", se = FALSE) +  # Add trend lines
    labs(title = "Accuracy Across Iterations for Each Model",
         x = "Iteration",
         y = "Accuracy") +
      theme(axis.text = element_text(size = 20, colour = "black"),
            axis.title.x = element_text(size = 22, face = "bold", margin = margin(t=10)),
            axis.title.y = element_text(size = 22, face = "bold", margin = margin(r=10)),
            legend.text = element_text(size = 18),
            legend.title = element_text(size = 20),
            panel.grid = element_blank(),
            panel.background = element_blank()))
  
  # ggplot2 Faceted Plot
  print(
    ggplot(performance_df %>% filter(Metric == "Accuracy"), aes(x = Iteration, y = Value)) +
    geom_line(color = "steelblue") +
    geom_point(alpha = 0.5) +
    facet_wrap(~ Model, ncol = 2) +
    labs(title = "Accuracy Across Iterations by Model",
         x = "Iteration",
         y = "Accuracy") +
    theme_minimal())
  
  # ggplot2 Error Bar Plot
  print(
    ggplot(summary_stats %>% filter(Metric == "Accuracy"), aes(x = reorder(Model, Mean), y = Mean)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), 
                  width = 0.2) +
    coord_flip() +  # Flip coordinates for better readability
    labs(title = "Mean Accuracy with Standard Deviation by Model",
         x = "Model",
         y = "Mean Accuracy") +
    theme_minimal())
  
  # Combine hyperparameter data
  if (length(hyperparameters_list) > 0) {
    hyperparameters_df <- bind_rows(hyperparameters_list)
    
    # Ensure numeric conversion for hyperparameter values
    hyperparameter_cols <- setdiff(names(hyperparameters_df), c("Iteration", "Model"))
    hyperparameters_df[hyperparameter_cols] <- lapply(hyperparameters_df[hyperparameter_cols], as.numeric)
    
    # Merge performance data with hyperparameter data
    merged_df <- merge(performance_df, hyperparameters_df, by = c("Iteration", "Model"), all.x = TRUE)
    
    # Reshape the merged data for plotting
    melted_df <- melt(merged_df, 
                      id.vars = c("Iteration", "Model", "Metric", "Value"), 
                      variable.name = "Hyperparameter", 
                      value.name = "Hyperparameter_Value",
                      na.rm = TRUE)  # Remove rows with NA in hyperparameters
    
    # Ensure Metric column exists before filtering
    if ("Metric" %in% colnames(melted_df)) {
      # Plot Hyperparameter Effect on Accuracy
      print(
        ggplot(melted_df %>% filter(Metric == "Accuracy"), aes(x = Hyperparameter_Value, y = Value)) +
          geom_point(alpha = 0.6, color = "lightblue") +
          geom_smooth(method = "loess", se = FALSE, color = "tomato") +
          facet_wrap(~ Hyperparameter, scales = "free_x") +
          labs(title = "Effect of Hyperparameters on Accuracy",
               x = "Hyperparameter Value",
               y = "Accuracy") +
          theme_minimal()
      )
    } else {
      warning("Metric column not found in melted_df. Skipping hyperparameter plot.")
    }
  } else {
    warning("No hyperparameter tuning data found. Skipping hyperparameter analysis.")
  }
  
  if (return_plot_data) {
    return(list(
      hyperparameter_data = melted_df,
      performance_df = performance_df,
      summary_stats = summary_stats
    ))
  }
}


##############################################################################
#                           Extract feature importance                       #
##############################################################################
extract_feature_importance <- function(ml_results, plot_feature_importance = TRUE) {
  library(caret)
  
  models_of_interest <- names(ML_bagging_results$model_predictions)
  
  iter_models <- ml_results$iter_models
  n_models <- length(iter_models)
  
  # Initialize a list to store feature importance for each model type
  feature_importance_list <- list()
  
  for (model_name in models_of_interest) {
    feature_importance_list[[model_name]] <- list()
    
    for (iter in 1:n_models) {
      model <- iter_models[[iter]][[model_name]]
      
      # Check if the model exists and varImp is applicable
      if (!is.null(model)) {
        imp <- try(varImp(model, scale = FALSE), silent = TRUE)
        
        if (class(imp) != "try-error") {
          # Extract the importance scores
          imp_df <- imp$importance
          
          # Some models might return multiple columns (e.g., for classification)
          # We'll take the overall importance by averaging if multiple columns exist
          if (ncol(imp_df) > 1) {
            imp_overall <- rowMeans(imp_df, na.rm = TRUE)
          } else {
            imp_overall <- imp_df
          }
          
          # Store in the list
          feature_importance_list[[model_name]][[iter]] <- imp_overall
        }
      }
    }
    
    # Convert the list of importance vectors to a data frame
    if (length(feature_importance_list[[model_name]]) > 0) {
      imp_matrix <- do.call(cbind, feature_importance_list[[model_name]])
      imp_df <- as.data.frame(imp_matrix)
      
      # Compute the average importance across iterations
      avg_imp <- rowMeans(imp_df, na.rm = TRUE)
      
      # Sort the features by average importance
      avg_imp_sorted <- sort(avg_imp, decreasing = TRUE)
      
      # Store back
      feature_importance_list[[model_name]] <- avg_imp_sorted
    } else {
      feature_importance_list[[model_name]] <- NULL
      warning(paste("No feature importance available for model:", model_name))
    }
  }
  
  if (isTRUE(plot_feature_importance)) {
    # Loop through each model in the list and plot
    for (model_name in names(feature_importance_list)) {
      # Skip if no feature importance available
      if (is.null(feature_importance_list[[model_name]])) next
      
      # Extract named vector of importances
      importance_vec <- feature_importance_list[[model_name]]
      
      # Convert to data frame
      imp_df <- data.frame(
        Feature    = names(importance_vec),
        Importance = as.numeric(importance_vec),
        stringsAsFactors = FALSE
      )
      
      # (Optional) Subset to the top 20 most important for clarity
      top_n <- 35
      imp_df <- imp_df[order(imp_df$Importance, decreasing = TRUE), ]
      imp_df <- head(imp_df, top_n)
      
      # Create bar plot
      p <- ggplot(imp_df, aes(x = reorder(Feature, Importance), y = Importance)) +
        geom_bar(stat = "identity", fill = "steelblue") +
        coord_flip() +
        theme_minimal(base_size = 14) +
        labs(title = paste("Feature Importance -", model_name),
             x = "Feature",
             y = "Importance") +
        theme(panel.grid = element_blank(),
              axis.text = element_text(colour = "black", size = 20),
              axis.title.x = element_text(size = 22, margin = margin(t=10)),
              axis.title.y = element_text(size = 22, margin = margin(r=10)))
      
      print(p)  # Print the plot
    }
  }
  
  return(feature_importance_list)
}


##############################################################################
#                           Model performance on test set                    #
##############################################################################
model_metrics <- function(feature_matrix, ML_bagging_results, plot_results = TRUE) {
  models <- names(ML_bagging_results$average_predictions)
  library(pROC)
  library(caret)
  
  # Append prediction scores to feature matrix
  feature_data_scores_appended <- feature_matrix
  for (model_name in models) {
    feature_data_scores_appended[, paste0("Prediction_Score_", model_name)] <- ML_bagging_results$average_predictions[[model_name]]
  }
  
  # Initialize a list to store confusion matrices and metrics
  model_metrics <- list()
  
  for (model_name in models) {
    
    # Extract the probability scores
    prob_scores <- feature_data_scores_appended[[paste0("Prediction_Score_", model_name)]]
    
    # Compute the ROC curve and AUC (which you already do)
    suppressMessages({
      roc_curve <- roc(response = feature_data_scores_appended$clinical, predictor = prob_scores)
      auc_value <- auc(roc_curve)
    })
    
    # Set threshold
    threshold <- 0.5
    
    # Convert probabilities to class predictions
    pred_class <- ifelse(prob_scores >= threshold, 1, 0)
    
    # Build a confusion matrix
    cm <- confusionMatrix(
      data = factor(pred_class, levels = c(0, 1)),  # predicted
      reference = factor(feature_data_scores_appended$clinical, levels = c(0, 1))  # actual
    )
    
    # Extract metrics from the confusion matrix
    accuracy <- cm$overall["Accuracy"]
    kappa <- cm$overall["Kappa"]
    sensitivity <- cm$byClass["Sensitivity"]     # same as recall
    specificity <- cm$byClass["Specificity"]
    precision <- cm$byClass["Pos Pred Value"]
    recall <- sensitivity                   # rename for clarity
    f1_score <- 2 * (precision * recall) / (precision + recall)
    balanced_acc <- cm$byClass["Balanced Accuracy"]
    
    # Store results
    model_metrics[[model_name]] <- list(
      ConfusionMatrix = cm,
      AUC = auc_value,
      Accuracy = accuracy,
      Kappa = kappa,
      Sensitivity = sensitivity,
      Specificity = specificity,
      Precision = precision,
      Recall = recall,
      F1_Score = f1_score,
      Balanced_Acc = balanced_acc
    )
  }
  
  if (isTRUE(plot_results)) {
    library(tidyverse)
    library(ggplot2)
    library(RColorBrewer)
    
    # Convert to a data frame, excluding the ConfusionMatrix object
    df <- map_dfr(
      names(model_metrics), 
      function(model_name) {
        x <- model_metrics[[model_name]]
        # Make a 1-row data frame for each model with numeric metrics
        data.frame(
          Model = model_name,
          AUC = as.numeric(x$AUC),
          Accuracy = as.numeric(x$Accuracy),
          Kappa = as.numeric(x$Kappa),
          Sensitivity = as.numeric(x$Sensitivity),
          Specificity = as.numeric(x$Specificity),
          Precision = as.numeric(x$Precision),
          Recall = as.numeric(x$Recall),
          F1_Score = as.numeric(x$F1_Score),
          Balanced_Acc = as.numeric(x$Balanced_Acc)
        )
      }
    )
    
    # Now pivot to long format
    df_long <- df %>%
      pivot_longer(
        cols = c(AUC, Accuracy, Kappa, Sensitivity, Specificity, Precision, Recall, F1_Score, Balanced_Acc),
        names_to = "Metric",
        values_to = "Value"
      )
    
    # Faceted Plot
    p <- ggplot(df_long, aes(x = Model, y = Value, fill = Model)) +
      geom_col(position = "dodge") +
      facet_wrap(~ Metric, scales = "free_y") +
      labs(
        title = "Model Performance (clinical) Comparison",
        x = "Model",
        y = "Metric Value"
      ) +
      scale_fill_brewer(palette = "Set2") +
      theme_bw(base_size = 14) +
      theme(
        legend.position = "none",
        strip.text = element_text(face = "bold")
      )
    
    print(p)
    
    return(list(model_metrics = model_metrics,
                feature_data_scores_appended = feature_data_scores_appended,
                plot_data = df))
  } else {
    return(list(model_metrics = model_metrics,
                feature_data_scores_appended = feature_data_scores_appended))
  }
}

