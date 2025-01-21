


ML_bagging <- function(n_models, feature_matrix, positive_set, negative_pool) {
  
  suppressMessages({
  library(caret)
  library(progress)
  library(foreach)
  library(doParallel)
  library(doSNOW)
  library(pROC)
  })
  
  cl <- makeCluster(detectCores() - 2)
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
    iter = seq_len(n_models),
    .combine = rbind,
    .packages = c("caret", "pROC"),
    .options.snow = opts
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
    # Select only numeric columns
    training_features_numeric <- training_features[, sapply(training_features, is.numeric)]
    training_features_numeric <- scale(training_features_numeric)
    
    # Train control with 10-fold cross-validation
    train_control <- trainControl(method = "cv", number = 10, 
                                  summaryFunction = twoClassSummary, 
                                  classProbs = TRUE, allowParallel = FALSE)
    
    n_features <- ncol(feature_matrix)-3 # n_features - Protein and label cols
    # Hyperparameter grids
    tune_grids <- list(
      glmnet = expand.grid(alpha = c(0, 0.5, 1), lambda = 10^seq(-4, 1, length = 10)),
      rf = expand.grid(mtry = seq(round(sqrt(n_features)) - 3, round(sqrt(n_features)) + 3, by = 1)),
      svmRadial = expand.grid(sigma = c(0.01, 0.1, 1, 10), C = c(0.1, 1, 10, 100)),
      knn = expand.grid(k = seq(round(sqrt(n_features))-3, round(sqrt(n_features))+3, by = 1)),  
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
    
    #Train all models on this iteration's data
    this_iter_models <- list()
    for (model_name in names(tune_grids)) {
      set.seed(iter)
      
      # XGboost may print depreciation message still
      if (model_name == "xgbTree") {
        Sys.setenv(XGBOOST_LOG_LEVEL = "0")
        suppressMessages(suppressWarnings(
          fit <- train(
            x = training_features_numeric, 
            y = training_labels, 
            method = model_name, 
            trControl = train_control, 
            tuneGrid = tune_grids[[model_name]],
            metric = "ROC", 
            verbose = FALSE
          )
        ))
      } else {
        suppressMessages(suppressWarnings(
          fit <- train(
            x = training_features_numeric, 
            y = training_labels, 
            method = model_name, 
            trControl = train_control, 
            tuneGrid = tune_grids[[model_name]],
            metric = "ROC"
          )
        ))
      }
      
      # Store the trained model
      this_iter_models[[model_name]] <- fit
    } # training end
    
    # Store results for each model in a small data.frame
    iteration_results <- data.frame()
    
    for (model_name in names(this_iter_models)) {
      model_probs <- predict(
        this_iter_models[[model_name]],
        newdata = feature_matrix[, sapply(feature_matrix, is.numeric)],
        type = "prob"
      )[,"Yes"]  # Probability that label is "Yes"
      
      # ROC using the 'clinical' column as the outcome (assuming 0/1)
      roc_curve <- suppressMessages(
        roc(
          response  = feature_matrix$clinical,
          predictor = model_probs
        )
      )
      auc_value <- pROC::auc(roc_curve)
      
      iteration_results <- rbind(
        iteration_results, 
        data.frame(
          iteration = iter,
          model     = model_name,
          AUC       = as.numeric(auc_value)
        )
      )
    }
    
    # Return iteration_results to be combined by foreach .combine='rbind'
    iteration_results
  } # do par end
  stopCluster(cl)
  
} # function end


library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
feature_matrix <- read.table("data/feature_matrix.txt", sep = "\t", header = T)
training_data <- data_sets_from_TTD(feature_matrix, ensembl)

ML_bagging_results <- ML_bagging(feature_matrix = training_data$feature_matrix,
                                 positive_set = training_data$positive_set, 
                                 negative_pool = training_data$negative_pool, 
                                 n_models = 2)


library(tidyverse)
library(ggplot2)

# Summarize
summary_by_model <- results_across_iterations %>%
  group_by(model) %>%
  summarise(
    meanAUC = mean(AUC),
    sdAUC   = sd(AUC),
    minAUC  = min(AUC),
    maxAUC  = max(AUC),
    nIters  = n()
  ) %>%
  arrange(desc(meanAUC))

summary_by_model

# Boxplot of AUC by model
ggplot(results_across_iterations, aes(x=model, y=AUC)) +
  geom_boxplot() +
  theme_bw() +
  ggtitle("AUC Distribution Over Multiple Negative-Pooling Iterations") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))





library(caret)
library(biomaRt)
library(pROC)

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
feature_matrix <- read.table("data/feature_matrix.txt", sep = "\t", header = T)

training_data <- data_sets_from_TTD(feature_matrix, ensembl)
n_features <- ncol(feature_matrix)-1

# Resample negatives for each iteration
set.seed(123)
negative_set <- training_data$negative_pool[sample(seq_len(nrow(training_data$negative_pool)), 
                                     nrow(training_data$positive_set), replace = FALSE), ]
training_set <- rbind(training_data$positive_set, negative_set)
training_labels <- training_set$approved
training_labels <- factor(
  training_labels,          # e.g., c(0, 1, 0, 1, ...)
  levels = c(0, 1),         # Original numeric levels
  labels = c("No", "Yes")   # New valid factor labels
)

# "Yes" is the second level --> 'Yes' is the positive class

training_features <- training_set[, !names(training_set) %in% c("approved", "Protein", "clinical")]
# Select only numeric columns
training_features_numeric <- training_features[, sapply(training_features, is.numeric)]
#training_features_numeric <- scale(training_features_numeric)

# Train control with 10-fold cross-validation
train_control <- trainControl(method = "cv", number = 10, summaryFunction = twoClassSummary, classProbs = TRUE)

# Hyperparameter grids
tune_grids <- list(
  glmnet = expand.grid(alpha = c(0, 0.5, 1), lambda = 10^seq(-4, 1, length = 10)),
  rf = expand.grid(mtry = seq(round(sqrt(n_features)) - 3, round(sqrt(n_features)) + 3, by = 1)),
  svmRadial = expand.grid(sigma = c(0.01, 0.1, 1, 10), C = c(0.1, 1, 10, 100)),
  knn = expand.grid(k = seq(round(sqrt(n_features))-3, round(sqrt(n_features))+3, by = 1)),  
  nb = expand.grid(fL = c(0, 0.5, 1, 2), usekernel = c(TRUE, FALSE), adjust = c(0.5, 1, 2)),
  nnet = expand.grid(size = c(5, 10, 15), decay = c(0.001, 0.01, 0.1, 1)),
  xgbTree = expand.grid(
    nrounds = c(50, 100, 150),       # Number of boosting iterations
    eta = c(0.01, 0.1, 1),         # Learning rate
    max_depth = c(3, 6, 9),          # Max depth of trees
    gamma = c(0, 1),              # Minimum loss reduction for split
    colsample_bytree = c(0.5, 0.8, 1), # Fraction of features per tree
    min_child_weight = c(1, 5),  # Minimum sum of instance weight for split
    subsample = c(1)       # Fraction of training samples per tree
  )
)


# Initialise an empty list to store model results
results <- list()

# Train models
for (model in names(tune_grids)) {
  set.seed(123)
  start <- Sys.time()
  cat("Training", model, "model...\n")
  
  if (model == "xgbTree") {
    Sys.setenv(XGBOOST_LOG_LEVEL = "0") 
    suppressMessages(suppressWarnings(
      fit <- train(
        x = training_features_numeric, 
        y = as.factor(training_labels), 
        method = model, 
        trControl = train_control, 
        tuneGrid = tune_grids[[model]],
        verbose = FALSE,
        metric = "ROC"
      )
    ))
  } else {
  suppressMessages(suppressWarnings(
    fit <- train(
    x = training_features_numeric, 
    y = as.factor(training_labels), 
    method = model, 
    trControl = train_control, 
    tuneGrid = tune_grids[[model]],
    metric = "ROC"
    )
  ))
  }

  cat("Elapsed time:", round(as.numeric(difftime(Sys.time(), start, units = "mins")), 2), "mins\n")
  
  results[[model]] <- fit
}


# Compare model performance
resamples <- resamples(results)
summary(resamples)

# Visualise performance comparison
bwplot(resamples, metric = "ROC")



all_predictions <- list()
for (model_name in names(results)) {
  cat("Predicting with:", model_name, "\n")
  all_predictions[[model_name]] <- predict(results[[model_name]],
                                           newdata = training_data$feature_matrix[, sapply(training_data$feature_matrix, is.numeric)],
                                           type = "prob"
                                           )[,2]

  # Evaluate performance on clincal targets
  roc_curve <- suppressMessages({
    roc(training_data$feature_matrix$clinical, as.numeric(all_predictions[[model_name]]) - 1)
  })
  print(auc(roc_curve))
}





