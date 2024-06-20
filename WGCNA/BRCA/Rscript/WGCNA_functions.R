# perform quality control
filter_low_expr <- function(tumour_matrix, control_matrix, sep = F) {
  if (sep == F) {
    # perform filtering on combined tumour and control matrix
    data <- merge(tumour_matrix, control_matrix, by = "row.names")
    data <- column_to_rownames(data, var = "Row.names")
    
    group <- factor(c(rep(1, length(colnames(tumour_matrix))), rep(2, length(colnames(control_matrix)))))
    counts_filt <- filterByExpr(data, group = group)
    
    print(table(counts_filt))
    counts_filt <- data[counts_filt, ]
    
    cat("Performing GSG check: \n")
    gsg <- goodSamplesGenes(t(counts_filt))
    
    geneCeck <- table(gsg$goodGenes)
    sampleCheck <- table(gsg$goodSamples)
    
    
    if (any(names(geneCeck) == "FALSE")) {
      cat("There are still bad genes in the dataset")
    } else if(any(names(geneCeck) == "FALSE")) {
      cat("There are still bad samples in the datase")
    }
    else {
      cat("GSG all clear, all samples ok")
    }
    
    return(counts_filt)
  }
  else if (sep == T) {
    # perform filtering fucntion on tumour and control matrix separately 
    cat("Filtering low expression genes from tumour group: \n")
    counts_filt <- filterByExpr(tumour_matrix)
    
    print(table(counts_filt))
    counts_filt <- tumour_matrix[counts_filt, ]
    
    cat("Performing GSG check: \n")
    gsg <- goodSamplesGenes(t(counts_filt))
    
    geneCeck <- table(gsg$goodGenes)
    sampleCheck <- table(gsg$goodSamples)
    
    if (any(names(geneCeck) == "FALSE")) {
      cat("There are still bad genes in the dataset")
    } else if(any(names(geneCeck) == "FALSE")) {
      cat("There are still bad samples in the datase")
    }
    else {
      cat("GSG all clear, all samples ok \n\n")
    }
    filtered_samples <- list(tumour = counts_filt)
    
    cat("Filtering low expression genes from control group: \n")
    counts_filt <- filterByExpr(control_matrix)
    
    print(table(counts_filt))
    counts_filt <- control_matrix[counts_filt, ]
    
    cat("Performing GSG check: \n")
    gsg <- goodSamplesGenes(t(counts_filt))
    
    geneCeck <- table(gsg$goodGenes)
    sampleCheck <- table(gsg$goodSamples)
    
    if (any(names(geneCeck) == "FALSE")) {
      cat("There are still bad genes in the dataset")
    } else if(any(names(geneCeck) == "FALSE")) {
      cat("There are still bad samples in the datase")
    }
    else {
      cat("GSG all clear, all samples ok")
    }
    filtered_samples[["control"]] <- counts_filt
  }
  return(filtered_samples)
}


# perform variance stabilizing transformation and transpose matrix for WGCNA
vst_norm <- function(counts_df) {
  matrix <- as.matrix(counts_df)
  matrix_norm <- varianceStabilizingTransformation(matrix)
  df_nrom <- as.data.frame(matrix_norm)
  transposed_norm_data <- t(df_nrom)
  
  return(transposed_norm_data)
}


# function to pick soft thresholding power with optional variable to plot only if already run function
pick_power <- function(WGCNA_data, network_type, sft_data) {
  if (missing(sft_data)) {
    start_time <- Sys.time()
    power <- c(c(1:10), seq(from = 12, to = 50, by = 2))
    
    # Call the network topology analysis function
    sft <- pickSoftThreshold(WGCNA_data,
                             powerVector = power,
                             networkType = network_type,
                             blockSize = 45000,
                             verbose = 2)
    
    sft_data <- sft$fitIndices
    
    # Visualise to pick power
    a1 <- ggplot(sft_data, aes(Power, SFT.R.sq, label = Power)) +
      geom_point() +
      geom_text(nudge_y = 0.1) +
      geom_hline(yintercept = 0.8, color = 'red') +
      labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
      theme_classic()
    
    a2 <- ggplot(sft_data, aes(Power, mean.k., label = Power)) +
      geom_point() +
      geom_text(nudge_y = 0.1) +
      labs(x = 'Power', y = 'Mean Connectivity') +
      theme_classic()
    
    grid.arrange(a1, a2, nrow = 2)
    
    sft_data
    time_elapsed <- Sys.time() - start_time
    cat("\n Time elapsed: ")
    print(time_elapsed)
    
    return(sft)
  }
  else {
    sft <- sft_data$fitIndices
    sft
    
    a1 <- ggplot(sft, aes(Power, SFT.R.sq, label = Power)) +
      geom_point() +
      geom_text(nudge_y = 0.1) +
      geom_hline(yintercept = 0.8, color = 'red') +
      labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
      theme_classic()
    
    a2 <- ggplot(sft, aes(Power, mean.k., label = Power)) +
      geom_point() +
      geom_text(nudge_y = 0.1) +
      labs(x = 'Power', y = 'Mean Connectivity') +
      theme_classic()
    
    grid.arrange(a1, a2, nrow = 2)
  }
}


# function to create seperate adj matrix from combined tumour and control counts matrix
sep_adj_matrix <- function(WGCNA_data, tumour_expr_df, control_expr_df, power) {
  start_time <- Sys.time()
  
  count_df <- as.data.frame(t(WGCNA_data))
  
  tumour_counts <- count_df[, colnames(count_df) %in% colnames(tumour_expr_df)]
  tumour_matrix <- t(tumour_counts)
  
  control_counts <- count_df[, colnames(count_df) %in% colnames(control_expr_df)]
  control_matrix <- t(control_counts)
  
  cat("Creating tumour adjacency... \n")
  tumour_adj <- adjacency(tumour_matrix, power = power, type = "unsigned")
  cat("Creating control adjacency... \n")
  control_adj <- adjacency(control_matrix, power = power, type = "unsigned")
  
  time_elapsed <- Sys.time() - start_time
  cat("Done, time elapsed: ")
  print(time_elapsed)
  adjcencies <- list(tumour = tumour_adj,
                     control = control_adj)
  
  return(adjcencies)
}


# identify modules
network_modules <- function(WGCNA_data, Power) {
  start_time <- Sys.time()
  bwnet <- blockwiseModules(WGCNA_data,
                            maxBlockSize = 45000,
                            TOMType = "signed",
                            networkType = "unsigned",
                            power = Power,
                            mergeCutHeight = 0.25,
                            numericLabels = FALSE,
                            randomSeed = 1234,
                            verbose = 3,
                            saveTOMs = FALSE)
  elapsed_time <- Sys.time() - start_time
  cat("Elapsed time: ")
  print(elapsed_time)
  
  
  # Plot the dendrogram
  plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                      c("unmerged", "merged"),
                      dendroLabels = FALSE,
                      addGuide = TRUE,
                      hang= 0.03,
                      guideHang = 0.05)
  
  return(bwnet)
}

