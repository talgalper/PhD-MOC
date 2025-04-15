### MOC pipeline Functions ###


plot_PCA <- function(expr_data, sample_info, output_plot_data = TRUE, circle_clust = FALSE, label_group = NULL, shape = NULL) {
  library(edgeR)
  library(ggplot2)
  library(ggrepel)
  library(ggalt)
  library(RColorBrewer)
  library(rlang)
  
  # Convert expression data frame into CPM-normalized + transposed matrix
  PCA_data <- cpm(as.matrix(expr_data), log = TRUE)
  PCA_data <- t(PCA_data)
  
  pca <- prcomp(PCA_data, scale. = TRUE, center = TRUE)
  pca_data <- as.data.frame(pca$x)
  
  pca_var <- pca$sdev^2
  pca_var_perc <- round(pca_var/sum(pca_var) * 100, digits = 2)
  
  # Merge sample mapping with PCA data
  pca_data <- merge(pca_data, sample_info, by.x = "row.names", by.y = 1)
  
  # Create a custom colour palette for stages
  groups <- unique(sample_info[ ,2])
  num_colors <- length(groups)
  colours <- brewer.pal(n = num_colors, name = "Set1")
  names(colours) <- groups
  
  # Set up the base aesthetic mapping depending on whether a shape variable is provided
  base_aes <- if (!is.null(shape)) {
    aes(PC1, PC2, color = Classification, shape = !!sym(shape))
  } else {
    aes(PC1, PC2, color = Classification)
  }
  
  # Build the ggplot object based on whether we want cluster encircling or not
  if (isTRUE(circle_clust)) {
    PCA_plot <- ggplot(pca_data, base_aes) +
      geom_point(size = 4) +
      geom_encircle(aes(group = Classification), s_shape = 0, expand = 0.05, color = "black") +
      scale_color_manual(values = colours) +
      theme_bw() +
      labs(x = paste0('PC1: ', pca_var_perc[1], ' %'),
           y = paste0('PC2: ', pca_var_perc[2], ' %')) +
      theme(
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 18, colour = "black"),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        panel.grid = element_blank()
      )
  } else {
    PCA_plot <- ggplot(pca_data, base_aes) +
      geom_point(size = 4) +
      scale_color_manual(values = colours) +
      theme_bw() +
      labs(x = paste0('PC1: ', pca_var_perc[1], ' %'),
           y = paste0('PC2: ', pca_var_perc[2], ' %')) +
      theme(
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 18, colour = "black"),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        panel.grid = element_blank()
      )
  }
  
  # Optionally add labels for a specific group if label_group is provided
  if (!is.null(label_group)) {
    PCA_plot <- PCA_plot + 
      geom_text_repel(
        data = subset(pca_data, Classification == label_group),
        aes(label = Row.names),
        size = 5,
        show.legend = FALSE
      )
  }
  
  # Print and optionally return the plot and plot data
  print(PCA_plot)
  
  if (output_plot_data) {
    return(list(PCA_plot = PCA_plot, plot_data = pca_data, pca_var_perc = pca_var_perc))
  }
}





# unpaired DE analysis
DE_analysis <- function(counts_matrix, sample_info) {
  
  data <- DGEList(counts = counts_matrix, group = sample_info$Classification)
  
  cat("Normalising library size \n")
  
  data <- normLibSizes(data)
  
  design <- model.matrix(~Classification, data = sample_info)
  
  cat("Estimating dispersion \n")
  
  data <- estimateDisp(data)
  
  cat("Conduct quasi-likelihood F-test \n")
  
  fit <- glmQLFit(data, design)
  qlf <- glmQLFTest(fit)
  toptags <- topTags(qlf, n = Inf)
  
  dif_exp <- decideTests(qlf, p = 0.05, adjust = "fdr", lfc = 1)
  print(summary(dif_exp))
  dif_exp_genes <- rownames(data$counts)[as.logical(dif_exp)]
  hits <- toptags$table[toptags$table$FDR < 0.1, ]
  colnames <- colnames(hits)
  hits$gene_id <- rownames(hits)
  hits <- hits[,c("gene_id", colnames)]
  dif_exp <- hits[dif_exp_genes, ]
  cat("Total DE:", nrow(dif_exp))
  
  return(list(input_data = counts_matrix,
              sample_info = sample_info,
              hits = hits,
              dif_exp = dif_exp,
              toptags = toptags,
              DGEList = data,
              fit = fit,
              qlf = qlf,
              design = design))
}


id_annot <- function(ensembl, data, col_id = 1, input_type, convert_to) {
  library(tidyverse)
  
  if (is.character(data)) {
    ensembl_annot <- getBM(
      attributes = c(input_type, convert_to), 
      filters = input_type, 
      values = data, 
      mart = ensembl
    )
    
    if (convert_to %in% "description") {
      ensembl_annot$description <- gsub("\\[.*?\\]", "", ensembl_annot$description)
    }
    
    # Create a matrix of NA's with rows equal to the length of 'data' and
    missing_df <- data.frame(
      matrix(NA, nrow = length(data), ncol = length(c(input_type, convert_to))),
      stringsAsFactors = FALSE
    )
    # Set the column names to match the ones used in getBM:
    colnames(missing_df) <- c(input_type, convert_to)
    # Fill in the input identifier column with the values in 'data'
    missing_df[, input_type] <- data
    
    # Combine the returned annotations with the missing rows data frame
    data_annot <- rbind(ensembl_annot, missing_df)
    data_annot <- data_annot[!duplicated(data_annot[,input_type]), ]
    
  } else {
    if (col_id == 0) {
      data <- rownames_to_column(data)
      col_id <- 1
    }
    
    ensembl_annot <- getBM(
      attributes = c(input_type, convert_to), 
      filters = input_type, 
      values = data[, col_id], 
      mart = ensembl
    )
    
    if (convert_to %in% "description") {
      ensembl_annot$description <- gsub("\\[.*?\\]", "", ensembl_annot$description)
    }    
    
    data_annot <- merge(
      ensembl_annot, 
      data, 
      by.x = input_type, 
      by.y = colnames(data)[col_id], 
      all.y = TRUE
    )
  }
  
  return(data_annot)
}

















