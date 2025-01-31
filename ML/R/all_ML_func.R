
ML_bagging <- function(n_models, feature_matrix, positive_set, negative_pool, big_grid = FALSE, 
                       models = c("glmnet", "rf", "svmRadial", "knn", "nb", "nnet", "xgbTree")) {
  
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
#                           Run ML bagging function                          #
##############################################################################
library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
feature_matrix <- read.table("data/feature_matrix.txt", sep = "\t", header = T)
training_data <- data_sets_from_TTD(feature_matrix, ensembl)

start <- Sys.time()
ML_bagging_results <- ML_bagging(feature_matrix = training_data$feature_matrix,
                                 positive_set = training_data$positive_set, 
                                 negative_pool = training_data$negative_pool, 
                                 n_models = 100)
print(Sys.time() - start)
rm(start)

load("~/OneDrive - RMIT University/PhD/large_git_files/ML/ML_bagging_100itr.RData")

##############################################################################
#                                 Plot results                               #
##############################################################################
plot_results <- function(return_plot_data = FALSE) {
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
plot_result <- plot_results(return_plot_data = TRUE)


##############################################################################
#                           Extract feature importance                       #
##############################################################################

extract_feature_importance <- function(ml_results, models_of_interest = c("glmnet", "rf", "svmRadial", "knn", "nb", "nnet", "xgbTree")) {
  library(caret)
  
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
  
  return(feature_importance_list)
}

feature_importance <- extract_feature_importance(ML_bagging_results)
