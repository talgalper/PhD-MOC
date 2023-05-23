library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(edgeR)
library(reshape2)

stages <- c("stage_I", "stage_II", "stage_III", "stage_IV")

for (i in seq_along(stages)) {
  stage <- stages[i]
  print(paste0("Running edgeR on ", stage, " samples"))
  
  
  gtex_data <- read.table("gtex/gene_reads_2017-06-05_v8_ovary.gct", skip = 2)
  
  # set col names of gtex data
  colnames(gtex_data) <- gtex_data[1, ]
  gtex_data <- gtex_data[-1, ]
  
  # rename cols
  colnames(gtex_data)[2] <- "gene_id"
  colnames(gtex_data)[3] <- "gene_name"
  gtex_data <- gtex_data[, -1]
  
  # remove number after dot point in ensembl id
  gtex_data$gene_id <- gsub("\\.[0-9]*$", "", gtex_data$gene_id)
  
  # get rid of gene name column
  gtex_data <- gtex_data[, -2]
  
  kylie_data <- read.csv(paste0("kylie/", stage, "_master_df.csv"))
  colnames(kylie_data)[1] <- "gene_id"
  
  # combine into single df
  df <- inner_join(kylie_data, gtex_data, by = c("gene_id"))
  rownames(df) <- df$gene_id
  df <- df[, -1]
  # convert to numeric because for some reason they are as chr
  df[] <- apply(df, 2, as.numeric)
  
  # subset genes that have a cpm of greater than one in more than 50% of the samples
  keepTheseGenes <- (rowSums(cpm(df) > 1) >= ncol(df)/2)
  print("Kept genes")
  print(summary(keepTheseGenes))
  
  #beforeFiltering_plot <- df %>% 
  #  cpm(log = TRUE) %>% 
  #  melt %>% 
  #  dplyr::filter(is.finite(value)) %>% 
  #  ggplot(aes(x = value, colour = Var2)) +
  #  geom_density() + 
  #  guides(colour = FALSE) +
  #  ggtitle("A. Before filtering", subtitle = paste0(nrow(df), " genes")) +
  #  labs(x = "logCPM", y = "Density")
  #
  #afterFiltering_plot <- df %>% 
  #  cpm(log = TRUE) %>% 
  #  magrittr::extract(keepTheseGenes,) %>%
  #  melt %>% 
  #  dplyr::filter(is.finite(value)) %>% 
  #  ggplot(aes(x = value, colour = Var2)) +
  #  geom_density() + 
  #  guides(colour = FALSE) +
  #  ggtitle("B. After filtering", subtitle = paste0(table(keepTheseGenes)[[2]], " genes"))+
  #  labs(x = "logCPM", y = "Density")
  #
  #beforeFiltering_plot
  #afterFiltering_plot
  
  
  
  df <- df[keepTheseGenes, ]
  
  
  
  # Select the columns that start with "unstranded"
  cols <- grep("^GAMuT", colnames(df))
  
  # Assign the unstranded columns to the cancer group
  cancer <- colnames(df)[cols]
  
  # Assign the rest of the columns to the healty group
  healthy <- setdiff(colnames(df), cancer)
  
  # Create the group variable
  group <- factor(c(rep("cancer", length(cancer)), rep("healthy", length(healthy))))
  
  
  data <- DGEList(counts = df, group = group)
  
  
  design <- model.matrix(~group)
  
  # Estimate a common negative binomial dispersion parameter for a DGE dataset with a general experimental design
  common <- estimateGLMCommonDisp(data, design, verbose = T)
  
  # Estimate the abundance-dispersion trend by Cox-Reid approximate profile likelihood.
  trend <- estimateGLMTrendedDisp(common, design)
  
  # Compute an empirical Bayes estimate of the negative binomial dispersion parameter for each tag, 
  # with expression levels specified by a log-linear model.
  tagwise <- estimateGLMTagwiseDisp(trend, design)
  
  # Plot the genewise biological coefficient of variation (BCV) against gene abundance (in log2 counts per million).
  plotBCV(tagwise, main = stage)
  ggsave(paste0("edgeR_results/kylie_data/", stage, "/", stage, "_BCV_plot.png"), device = "png")
  
  # Fit a negative binomial generalized log-linear model to the read counts for each gene. 
  # Conduct genewise statistical tests for a given coefficient or coefficient contrast.
  fit <- glmFit(tagwise, design)
  
  # Fit a negative binomial generalized log-linear model to the read counts for each gene. 
  # Conduct genewise statistical tests for a given coefficient or coefficient contrast.
  lrt <- glmLRT(fit, coef = 2)
  
  # Extract the most differentially expressed genes (or sequence tags) from a test object, 
  # ranked either by p-value or by absolute log-fold-change.
  toptags <- topTags(lrt, n = Inf)
  
  # Identify which genes are significantly differentially expressed from 
  # an edgeR fit object containing p-values and test statistics.
  dif_exp <- decideTestsDGE(lrt, p = 0.05, adjust = "fdr", lfc = 1)
  print(summary(dif_exp))
  
  dif_exp_genes <- rownames(tagwise)[as.logical(dif_exp)]
  
  # Make a mean-difference plot of two libraries of count data with smearing of points with very low counts, 
  # especially those that are zero for one of the columns.
  plotSmear(lrt, de.tags = dif_exp_genes, main = stage)
  ggsave(paste0("edgeR_results/kylie_data/", stage, "/",  stage, "_smear_plot.png"), device = "png")
  
  # create a results df
  hits <- toptags$table[toptags$table$FDR < 0.1, ]
  colnames <- colnames(hits)
  hits$gene_id <- rownames(hits)
  hits <- hits[,c("gene_id", colnames)]
  
  # plot Pvalues of different logFC scores
  ggplot(hits, aes(x=logFC, y=PValue)) + geom_point() + labs(title = stage)
  ggsave(paste0("edgeR_results/kylie_data/", stage, "/", stage, "_Pvalue_plot.png"), device = "png")
  
  
  write.csv(hits, paste0("edgeR_results/kylie_data/", stage, "/", stage, "_edgeR_hits.csv"))
}


