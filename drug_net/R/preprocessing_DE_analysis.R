


#' combines TCGA and GTEx data for preprocessing
#' 
#' @param TCGA_data TCGA unstranded data. Output by `get_TCGA_RNAseq_data`
#' @param GTEx_data Data downloaded from the GTEx portal
#' @examples
#' GTEx_data <- read.table("bulk-gex_v8_rna-seq_counts-by-tissue_gene_reads_2017-06-05_v8_breast_mammary_tissue.gct", skip = 2, header = T)
#' data <- TCGA_GTEx_combine(TCGA_data, GTEx_data)
#' @return Returns a dataframe with the two input datas combined, ready for differential expression analysis
#' @export
TCGA_GTEx_combine <- function(TCGA_data, GTEx_data) {
  GTEx_raw <- GTEx_data
  GTEx_data <- GTEx_raw[, -c(1:3)]
  rownames(GTEx_data) <- GTEx_raw$Name
  rownames(GTEx_data) <- gsub("\\.\\d+", "", rownames(GTEx_data))
  
  data <- merge(TCGA_data, GTEx_data, by = "row.names")
  data <- column_to_rownames(data, "Row.names")
  
  return(data)
}




#' Removes low activity transcripts from data. Requires a file containing all samples including healthy/control and disease groups.
#' 
#' @param data data contianing healthy/control and disease groups
#' @param top_x_samples percentile range that a trasncript must appear in
#' @param min_sample Percentage value indicating the number of samples the trasncript must appear within the range specified in `top_x_samples`
#' @examples
#' treated_data <- remove_low_activity_genes(data, top_x_samples = 0.50, min_samples = 0.10)
#' @return data frame with low activity genes remove.
#' @export
remove_low_activity_genes <- function(data, top_x_samples, min_samples) {
  
  num_genes <- nrow(data)
  top_threshold <- ceiling(num_genes * 0.50)
  
  # min required samples
  num_samples <- ncol(data)
  min_required_samples <- ceiling(num_samples * 0.10)
  
  # Extract the gene names and data columns
  gene_names <- rownames(data)
  data_columns <- data[, 2:ncol(data)]
  
  # Initialise an empty data frame
  results_df <- data.frame(gene_id = gene_names)
  
  # Iterate through each data column and perform the test
  for (col_id in seq_along(data_columns)) {
    # Get the column name and the data
    col_name <- colnames(data_columns)[col_id]
    col_data <- data_columns[, col_id]
    
    # Sort the data and mark genes within the top_threshold rows as TRUE
    sorted_indices <- order(col_data, decreasing = TRUE)
    top_indices <- sorted_indices[1:top_threshold]
    
    # Create a logical vector for gene presence in top_threshold
    gene_presence <- rep(FALSE, nrow(data))
    gene_presence[top_indices] <- TRUE
    
    # Add to the results data frame
    results_df <- cbind(results_df, gene_presence)
  }
  
  colnames(results_df) <- c("gene_id", colnames(data_columns))
  
  # Calculate the number of TRUE values for each gene across columns
  gene_counts <- rowSums(results_df[, -1]) # Exclude the first column ("gene_id")
  failed_genes <- gene_counts < min_required_samples
  cat("Failed genes", "\n")
  print(summary(failed_genes))
  
  results <- subset(data, !(rownames(data) %in% gene_names[failed_genes]))
  
  return(results)
}


#' Performs differential expression analysis on TCGA data
#' 
#' @param data Raw counts from TCGA and healthy cell group
#' @param show_plots Displays BCV, smear and Pvalu plots
#' @examples
#' # example code
#' DE_results <- TCGA_DE_analysis(data)
#' @return Dataframe of significantly differentially expressed genes as well as plots.
#' @export
TCGA_DE_analysis <- function(data, show_plots = TRUE) {
  
  cat("Make take a moment...")
  
  # Select the columns that start with "TCGA"
  cols <- grep("^TCGA", colnames(data))
  
  # Assign the unstranded columns to the cancer group
  cancer <- colnames(data)[cols]
  
  # Assign the rest of the columns to the healty group
  benign <- setdiff(colnames(data), cancer)
  
  # Create the group variable
  group <- factor(c(rep("cancer", length(cancer)), rep("benign", length(benign))))
  
  data <- DGEList(counts = data, group = group)
  
  design <- model.matrix(~group)
  
  # Estimate a common negative binomial dispersion parameter for a DGE dataset with a general experimental design
  common <- estimateGLMCommonDisp(data, design, verbose = T)
  
  # Estimate the abundance-dispersion trend by Cox-Reid approximate profile likelihood.
  trend <- estimateGLMTrendedDisp(common, design)
  
  # Compute an empirical Bayes estimate of the negative binomial dispersion parameter for each tag, 
  # with expression levels specified by a log-linear model.
  tagwise <- estimateGLMTagwiseDisp(trend, design)
  
  # Fit a negative binomial generalized log-linear model to the read counts for each gene. 
  # Conduct genewise statistical tests for a given coefficient or coefficient contrast.
  fit <- glmFit(tagwise, design)
  
  # Conduct genewise statistical tests for a given coefficient or coefficient contrast.
  lrt <- glmLRT(fit, coef = 2)
  
  # Extract the most differentially expressed genes (or sequence tags) from a test object, 
  # ranked either by p-value or by absolute log-fold-change.
  toptags <- topTags(lrt, n = Inf)
  
  # Identify which genes are significantly differentially expressed from 
  # an edgeR fit object containing p-values and test statistics.
  dif_exp <- decideTestsDGE(lrt, p = 0.05, adjust = "fdr", lfc = 2)
  print(summary(dif_exp))
  
  dif_exp_genes <- rownames(tagwise)[as.logical(dif_exp)]
  
  # create a results df
  hits <- toptags$table[toptags$table$FDR < 0.1, ]
  colnames <- colnames(hits)
  hits$gene_id <- rownames(hits)
  hits <- hits[,c("gene_id", colnames)]
  
  dif_exp <- hits[dif_exp_genes, ]
  
  if (show_plots == TRUE) {
    # Plot the genewise biological coefficient of variation (BCV) against gene abundance (in log2 counts per million).
    plotBCV(tagwise, main = "biological coefficient of variation")
    
    # Make a mean-difference plot of two libraries of count data with smearing of points with very low counts, 
    # especially those that are zero for one of the columns.
    plotSmear(lrt, de.tags = dif_exp_genes, main = "Mean-difference plot")
    
    # plot Pvalues of different logFC scores
    ggplot(hits, aes(x=logFC, y=-log(FDR))) + geom_point() + labs(title = "Adjusted logFC")
  }
  
  return(list(hits, dif_exp))
}


#' Convert IDs using biomaRt human ensembl database
#' 
#' @param gene_list A vector containing the terms to be converted
#' @param from The format of the input data
#' @param to The format of the desired output
#' @param description If changed to TRUE will include an addtional column with description of term
#' @examples
#' ID_converted <- convert_ID(gene_list, from = "ensembl_gene_id", to = "external_gene_name")
#' @return description
#' @export
convert_ID <- function(gene_list, from, to, description = FALSE) {
  
  if (!exists("ensembl")) {
    ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  }
  
  attributes <- listAttributes(ensembl)
  
  if (from %in% attributes$name) {
    
    id_convert <- getBM(attributes = c(from, to), 
                        filters = from, 
                        values = gene_list, 
                        mart = ensembl)
    
    id_convert <- subset(id_convert, to != "")
    
    # check for duplicate uniprot ids
    id_convert <- distinct(id_convert)
    
    missing_terms <- setdiff(gene_list, id_convert[[from]])
    
    if (description == TRUE) {
      # add gene description
      gene_description <- getBM(attributes = c(to, "description"), 
                                filters = to, 
                                values = id_convert[[to]], 
                                mart = ensembl)
      
      id_description <- merge(id_convert, gene_description, by = to)
      id_description$description <- gsub("\\s*\\[.*?\\]", "", id_description$description)
      
      if (length(missing_terms != 0)){
        cat("Could not be converted:", paste(missing_terms, collapse = ", "))
      }
      
      return(id_description)
    }
    
    if (length(missing_terms != 0)){
      cat("Could not be converted:", paste(missing_terms, collapse = ", "))
    }
    
    return(id_convert)
    
  } else {
    cat("ID name not found. Run following code for attribute list: ", "\n")
    cat("ensembl <- useMart('ensembl', dataset = 'hsapiens_gene_ensembl')", "\n")
    cat("attributes <- listAttributes(ensembl)")
  }
}