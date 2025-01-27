
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
  models <- c("glmnet", "rf", "svmRadial", "knn", 
              # "nb", 
              "nnet"#, 
              # "xgbTree"
              )
  model_predictions <- list()
  for (model in models) {
    model_predictions[[model]] <- matrix(0, nrow = nrow(feature_matrix), ncol = n_models)
  }

  # Set up parallel backend using doSNOW and 80% of cpu capacity
  
  cl <- makeSOCKcluster(round(detectCores()*0.8))
  registerDoSNOW(cl)
  on.exit(stopCluster(cl))
  
  
  # # Define the silencing function 
  # silence_all_output_null <- function(expr) {
  #   # Determine the null device based on the operating system
  #   nullfile <- if (.Platform$OS.type == "windows") "NUL" else "/dev/null"
  #   
  #   # Open connections to the null device
  #   conn_out <- file(nullfile, open = "w")
  #   conn_msg <- file(nullfile, open = "w")
  #   
  #   # Divert standard output and messages to the null device
  #   sink(conn_out)
  #   sink(conn_msg, type = "message")
  #   
  #   # Ensure that sinks are always restored, even if an error occurs
  #   on.exit({
  #     # Restore message sink first
  #     if (sink.number(type = "message") > 0) sink(type = "message")
  #     # Restore standard output sink
  #     if (sink.number() > 0) sink()
  #     # Close the connections
  #     close(conn_msg)
  #     close(conn_out)
  #   }, add = TRUE)
  #   
  #   # Muffle warnings within the silenced block
  #   withCallingHandlers(
  #     expr,
  #     warning = function(w) invokeRestart("muffleWarning")
  #   )
  #   
  #   invisible(NULL)
  # }
  
  # Export necessary variables and functions to the workers
  clusterExport(cl, 
                list = c("numeric_features", "feature_matrix", "positive_set", "negative_pool", "models"), 
                envir = environment())
  
  # setup progress bar for parallelisation
  pb <- txtProgressBar(max = n_models, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  parallel_result <- foreach(
    iter = 1:n_models, 
    #.options.snow = opts, 
    .combine = "list",
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
      #nb = expand.grid(fL = c(0, 0.5, 1, 2), usekernel = c(TRUE, FALSE), adjust = c(0.5, 1, 2)),
      nnet = expand.grid(size = c(5, 10, 15), decay = c(0.001, 0.01, 0.1, 1))#,
      # xgbTree = expand.grid(
      #   nrounds = c(50, 100, 150),   
      #   eta = c(0.01, 0.1, 1),    
      #   max_depth = c(3, 6, 9),   
      #   gamma = c(0, 1),      
      #   colsample_bytree = c(0.5, 0.8, 1), 
      #   min_child_weight = c(1, 5),  
      #   subsample = c(1)       
      # )
    )
    
    # initialise local data
    local_model_predictions <- list()
    # Train all models on this iteration's data within the silenced block
    this_iter_models <- list()
    #silence_all_output_null({
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
            verbose = FALSE
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
    #})
    
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
    list(iter_models = this_iter_models,
         model_predictions = local_model_predictions)
  } # dopar end
  
  stopCluster(cl)
  close(pb)
  
  iter_models <- list()
  # Combine the foreach results into model_predictions and resamples
  for (iter in 1:n_models) {
    iteration_result <- parallel_result[[iter]]
    iter_models[[iter]] <- parallel_result[[iter]][["iter_models"]]
    
    # Assign predictions
    for (model_name in models) {
      model_predictions[[model_name]][, iter] <- iteration_result$model_predictions[[model_name]]
    }
  }
  
  # mean prediction across iterations
  average_predictions <- lapply(model_predictions, function(pred_matrix) {
    rowMeans(pred_matrix)
  })
  
  return(list(model_predictions = model_predictions,
              average_predictions = average_predictions))
} # function end



library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
feature_matrix <- read.table("data/feature_matrix.txt", sep = "\t", header = T)
training_data <- data_sets_from_TTD(feature_matrix, ensembl)

# Select only numeric columns
#training_data$feature_matrix_numeric <- training_data$feature_matrix[, sapply(training_data$feature_matrix[,c(2:70)], is.numeric)]

start <- Sys.time()
ML_bagging_results <- ML_bagging(feature_matrix = training_data$feature_matrix,
                                 positive_set = training_data$positive_set, 
                                 negative_pool = training_data$negative_pool, 
                                 n_models = 2)
print(Sys.time() - start)
rm(start)
