
ML_bagging <- function(n_models, feature_matrix, positive_set, negative_pool) {
  
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
  
  # Initialise matrix for each model
  models <- c("glmnet", "rf", "svmRadial", "knn", "nb", "nnet", "xgbTree")
  
  model_predictions <- list()
  for (model in models) {
    model_predictions[[model]] <- matrix(0, nrow = nrow(feature_matrix), ncol = n_models)
  }
  
  # Initialize resamples list
  resamples_list <- vector("list", n_models)
  
  # Set up parallel backend using doSNOW and 80% of cpu capacity
  cl <- makeSOCKcluster(min(round(detectCores()*0.5), n_models))
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
                                  summaryFunction = twoClassSummary, 
                                  classProbs = TRUE, allowParallel = FALSE)
    
    n_features <- ncol(feature_matrix) - 3 # Adjust as per your data
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
    
    # initialise local data
    local_model_predictions <- list()
    # Train all models on this iteration's data within the silenced block
    this_iter_models <- list()
    for (model_name in names(tune_grids)) {
      set.seed(iter)
      
      if (model_name == "xgbTree") {
        Sys.setenv(XGBOOST_LOG_LEVEL = "0")
        fit <- train(
          x = training_features_numeric, 
          y = training_labels, 
          method = model_name, 
          trControl = train_control, 
          tuneGrid = tune_grids[[model_name]],
          metric = "ROC", 
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
          metric = "ROC"
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


##############################################################################
#                                 Plot results                               #
##############################################################################
plot_results <- function(return_plot_data = FALSE) {
  library(ggplot2)
  
  # Initialize a data frame to store ROC for each model and iteration
  performance_df <- data.frame(
    Iteration = integer(),
    Model = character(),
    ROC = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Loop through each iteration and extract ROC
  for (i in 1:length(ML_bagging_results$resamples)) {
    resample <- ML_bagging_results$resamples[[i]]
    # Extract ROC for each model
    roc_values <- summary(resample)$statistics$ROC[, "Mean"]
    
    # Create a temporary data frame for this iteration
    temp_df <- data.frame(
      Iteration = i,
      Model = names(roc_values),
      ROC = as.numeric(roc_values),
      stringsAsFactors = FALSE
    )
    
    # Append to the main data frame
    performance_df <- rbind(performance_df, temp_df)
  }
  
  # Compute summary statistics
  summary_stats <- aggregate(ROC ~ Model, data = performance_df, 
                             FUN = function(x) c(Mean = mean(x), SD = sd(x)))
  
  # Clean up the summary_stats data frame
  summary_stats <- do.call(data.frame, summary_stats)
  names(summary_stats) <- c("Model", "Mean_ROC", "SD_ROC")
  
  # Base R Boxplot
  boxplot(ROC ~ Model, data = performance_df,
          main = "Distribution of ROC AUC Across Iterations",
          ylab = "ROC AUC",
          col = "lightblue",
          las = 2)  # Rotate x-axis labels for better readability
  
  # ggplot2 Line Plot
  ggplot(performance_df, aes(x = Iteration, y = ROC, color = Model)) +
    geom_line(alpha = 0.6) +
    geom_smooth(method = "loess", se = FALSE) +  # Add trend lines
    labs(title = "ROC AUC Across Iterations for Each Model",
         x = "Iteration",
         y = "ROC AUC") +
    theme_minimal()
  
  # ggplot2 Faceted Plot
  ggplot(performance_df, aes(x = Iteration, y = ROC)) +
    geom_line(color = "blue") +
    geom_point(alpha = 0.5) +
    facet_wrap(~ Model, ncol = 2) +
    labs(title = "ROC AUC Across Iterations by Model",
         x = "Iteration",
         y = "ROC AUC") +
    theme_minimal()
  
  
  # Convert average_predictions to a data frame
  average_predictions_df <- data.frame(
    Model = names(ML_bagging_results$average_predictions),
    Average_Prediction = unlist(ML_bagging_results$average_predictions),
    stringsAsFactors = FALSE
  )
  
  # ggplot2 Error Bar Plot
  ggplot(summary_stats, aes(x = reorder(Model, Mean_ROC), y = Mean_ROC)) +
    geom_bar(stat = "identity", fill = "lightblue") +
    geom_errorbar(aes(ymin = Mean_ROC - SD_ROC, ymax = Mean_ROC + SD_ROC), 
                  width = 0.2) +
    coord_flip() +  # Flip coordinates for better readability
    labs(title = "Mean ROC AUC with Standard Deviation by Model",
         x = "Model",
         y = "Mean ROC AUC") +
    theme_minimal()
  
  
  # analyse hyperparamters
  library(reshape2)
  # Initialize a list to store hyperparameters
  hyperparameters_list <- list()
  
  for (i in 1:length(ML_bagging_results$iter_models)) {
    iter_models <- ML_bagging_results$iter_models[[i]]
    for (model_name in names(iter_models)) {
      model <- iter_models[[model_name]]
      best_tune <- model$bestTune
      # Add Iteration and Model columns
      best_tune$Iteration <- i
      best_tune$Model <- model_name
      # Append to the list
      hyperparameters_list[[length(hyperparameters_list) + 1]] <- best_tune
    }
  }
  
  # Identify all unique hyperparameter names across all models
  all_hyperparams <- unique(unlist(lapply(hyperparameters_list, names)))
  # Remove 'Iteration' and 'Model' from hyperparameter names
  all_hyperparams <- setdiff(all_hyperparams, c("Iteration", "Model"))
  
  # Function to standardize hyperparameter data frames
  standardize_hyperparams <- function(df, all_params) {
    missing_params <- setdiff(all_params, names(df))
    if(length(missing_params) > 0){
      df[missing_params] <- NA
    }
    # Ensure the order of columns is consistent
    df <- df[, c("Iteration", "Model", all_params)]
    return(df)
  }
  
  # Apply the standardization to each hyperparameter data frame
  hyperparameters_list_standardized <- lapply(hyperparameters_list, standardize_hyperparams, all_params = all_hyperparams)
  
  # Combine the standardized data frames into one
  hyperparameters_df <- do.call(rbind, hyperparameters_list_standardized)
  
  # Convert appropriate columns to numeric (if they aren't already)
  for(param in all_hyperparams){
    hyperparameters_df[[param]] <- as.numeric(as.character(hyperparameters_df[[param]]))
  }
  
  merged_df <- merge(performance_df, hyperparameters_df, by = c("Iteration", "Model"))
  
  # Melt the hyperparameters (excluding Iteration, Model, ROC)
  melted_df <- melt(merged_df, 
                    id.vars = c("Iteration", "Model", "ROC"), 
                    variable.name = "Hyperparameter", 
                    value.name = "Value",
                    na.rm = TRUE)  # Remove rows with NA in hyperparameters
  
  # Convert Value to numeric (if not already)
  melted_df$Value <- as.numeric(as.character(melted_df$Value))
  
  # Create a combined factor for Model and Hyperparameter for faceting
  melted_df$Model_Hyperparameter <- paste(melted_df$Model, melted_df$Hyperparameter, sep = "_")
  # Plot
  ggplot(melted_df, aes(x = Value, y = ROC)) +
    geom_point(alpha = 0.6, color = "lightblue") +
    geom_smooth(method = "loess", se = FALSE, color = "tomato") +
    facet_wrap(~ Model_Hyperparameter, scales = "free_x") +
    labs(title = "Effect of Hyperparameters on ROC AUC",
         x = "Hyperparameter Value",
         y = "ROC AUC") +
    theme_minimal()
  
  if (isTRUE(return_plot_data)) {
    return(list(hyperparameter_data = melted_df,
                performance_df = performance_df,
                summary_stats = summary_stats))
  }
}




