### edgeR wrapper functions for DE ###


# unpaired quality control
filter_low_expr <- function(disease_data, control_data) {
  merged_df <- merge(control_data, disease_data, by = "row.names")
  merged_df <- column_to_rownames(merged_df, "Row.names")
  
  sample_info <- data.frame(
    sample = colnames(merged_df),
    group = c(rep("control", ncol(control_data)), rep("disease", ncol(disease_data))))
  sample_info$group <- factor(sample_info$group, levels = c("control", "disease"))
  
  counts_filt <- filterByExpr(merged_df, group = sample_info$group)
  counts_filt <- merged_df[counts_filt, ]
  
  low_exp_genes <- merged_df[!rownames(merged_df) %in% rownames(counts_filt), ]
  
  return(list(counts_filt = counts_filt,
              sample_info = sample_info,
              low_exp_genes = low_exp_genes))
}

# paired quality control
paired_filter_low_expr <- function(disease_data, control_data, paired_info) {
  subtype <- colnames(paired_info)
  subtype <- subtype[3]
  
  merged_df <- merge(control_data, disease_data, by = "row.names")
  merged_df <- column_to_rownames(merged_df, "Row.names")
  
  merged_df <- merged_df[colnames(merged_df) %in% c(paired_info$normal, paired_info[[subtype]])] # subset counts for paired samples
  
  sample_info <- data.frame(samples = c(paired_info$normal, paired_info[[3]]),
                            pat_id = c(paired_info$bcr_patient_barcode, paired_info$bcr_patient_barcode),
                            group = c(rep("control", length(paired_info$normal)), rep("disease", length(paired_info[[3]]))))
  sample_info$group <- factor(sample_info$group, levels = c("control", "disease"))
  sample_info$pat_id <- factor(sample_info$pat_id)
  
  counts_filt <- filterByExpr(merged_df, group = sample_info$group)
  counts_filt <- merged_df[counts_filt, ]
  
  low_exp_genes <- merged_df[!rownames(merged_df) %in% rownames(counts_filt), ]
  
  return(list(counts_filt = counts_filt,
              sample_info = sample_info,
              low_exp_genes = low_exp_genes))
}

# unpaired DE analysis
DE_analysis <- function(counts_matrix, sample_info) {
  
  data <- DGEList(counts = counts_matrix, group = sample_info$group)
  
  cat("Normalising library size \n")
  
  data <- normLibSizes(data)
  
  design <- model.matrix(~group, data = sample_info)
  
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
  
  return(list(hits = hits,
              dif_exp = dif_exp,
              toptags = toptags,
              data = data,
              fit = fit,
              qlf = qlf,
              design = design))
}

# paired DE analysis
paired_DE_analysis <- function(counts_matrix, sample_info) {
  
  data <- DGEList(counts = counts_matrix, group = sample_info$group)
  
  cat("Normalising library size \n")
  
  data <- normLibSizes(data)
  
  design <- model.matrix(~sample_info$pat_id+sample_info$group)
  
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
  
  return(list(hits = hits,
              dif_exp = dif_exp,
              toptags = toptags,
              data = data,
              fit = fit,
              qlf = qlf,
              design = design))
}


# print DE counts
print_summary <- function(DE_results) {
  qlf <- DE_results$qlf
  data <- DE_results$data
  
  dif_exp <- decideTests(qlf, p = 0.05, adjust = "fdr", lfc = 1)
  print(summary(dif_exp))
  
  toptags <- topTags(qlf, n = Inf)
  dif_exp_genes <- rownames(data$counts)[as.logical(dif_exp)]
  hits <- toptags$table[toptags$table$FDR < 0.1, ]
  colnames <- colnames(hits)
  hits$gene_id <- rownames(hits)
  hits <- hits[,c("gene_id", colnames)]
  dif_exp <- hits[dif_exp_genes, ]
  cat("total:", nrow(dif_exp))
}

# returns a warning for some reason
plot_DE_results <- function(DE_results) {
  # Plot the genewise biological coefficient of variation (BCV) against gene abundance (in log2 counts per million).
  plotBCV(DE_results$data, main = "biological coefficient of variation")
  
  # Make a mean-difference plot of two libraries of count data with smearing of points with very low counts, 
  # especially those that are zero for one of the columns.
  plotSmear(DE_results$qlf, de.tags = DE_results$dif_exp_genes, main = "Mean-difference plot")
  
  # plot Pvalues of different logFC scores
  ggplot(DE_results$hits, aes(x=logFC, y=-log(FDR))) + geom_point() + labs(title = "Adjusted logFC")
}
