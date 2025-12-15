### Updated analysis of individual BRCA subtypes ###
### This analysis is in response to thesis reveiwer feedback ###

args <- commandArgs(trailingOnly = TRUE)


## function to analyse each subtype
subtype_analysis <- function(subtype_data, control_data, subtype_name, PCA = TRUE, run_PCSF = FALSE) {
  suppressMessages(library(edgeR))
  # laod in good DE functions from MOC_pipe
  source("../MOC_pipe/R/functions.R")
  
  # create sample info
  sample_info <- data.frame(sample = c(colnames(subtype_data), colnames(control_data)),
                            group = c(rep(subtype_name, ncol(subtype_data)), rep("control", ncol(control_data))))
  
  # combine sample type data
  load("RData/TCGA_query.RData")
  common <- common[, c(1,3)]
  sample_info <- merge(sample_info, common, by.x = "sample", by.y = "cases", all.x = T)
  sample_info$sample_type <- ifelse(is.na(sample_info$sample_type), "Healthy", sample_info$sample_type)
  
  # check for and remove STN samples
  STN_samples <- sample_info$sample[sample_info$sample_type == "Solid Tissue Normal"]
  if (length(STN_samples) > 0) {
    cat("Removing STN samples:", length(STN_samples), "\n")
    subtype_data <- subtype_data[, !colnames(subtype_data) %in% STN_samples]
    sample_info <- sample_info[!sample_info$sample %in% STN_samples, ] # remove STN samples from sample_info
  }
  
  # create DE dataset
  merged_df <- merge(control_data, subtype_data, by = "row.names")
  merged_df <- tibble::column_to_rownames(merged_df, "Row.names")
  
  # filter low counts
  counts_filt <- filterByExpr(merged_df, group = sample_info$group)
  print(table(counts_filt))
  counts_filt <- merged_df[counts_filt, ]
  low_exp_genes <- merged_df[!rownames(merged_df) %in% rownames(counts_filt), ]
  
  # optional PCA data visualisation
  if (PCA == TRUE) {
    cat("Plotting PCA...", "\n")
    PCA_data <- plot_PCA(expr_data = counts_filt, 
                         sample_info = sample_info, 
                         colour = "group", 
                         shape = "sample_type")
  } else {
    PCA_data <- NULL
  }
  
  # offer abort if PCA doest look good
  # ans <- readline(prompt = "Does everything look ok? continue... (y/n): ")
  # if (tolower(ans) != "y") {
  #   cat("Aborting.\n")
  #   return(invisible(NULL))
  # }
  
  # DE analysis
  cat("Performing DE analysis...", "\n")
  DE_data <- DE_analysis(counts_matrix = counts_filt,
                         sample_info = sample_info, 
                         group = "group")
  

  # optionally run PCSF using DE scores
  if (run_PCSF == TRUE) {
    cat("\n", "Running PCSF...", "\n")
    suppressMessages(library(PCSF))
    # load in STRING data
    STRING_edge <- read.csv("latest_run/intermediate/STRING_physical_ENSG.csv")
    STRING_edge <- STRING_edge[!duplicated(t(apply(STRING_edge, 1, sort))), ]
    # set seed for reproducibility 
    set.seed(1234)
    # construct interactome
    ppi <- construct_interactome(STRING_edge)
    # set terminals
    dif_exp <- DE_data$dif_exp
    dif_exp$logFC_abs <- abs(dif_exp$logFC)
    terminals <- setNames(as.numeric(dif_exp$logFC), dif_exp$gene_id)
    # run PCSF with random noise
    start_time <- Sys.time()
    subnet <- PCSF_rand(ppi, terminals, n = 50, r = 0.1, b = 1, w = 2, mu = 0.0005)
    elapsed_time <- Sys.time() - start_time
    print(elapsed_time)
    # plot subnetwork
    plot.PCSF(subnet, node_label_cex = 15)
    
    # extract cluster data
    clust <- components(subnet)
    df <- data.frame(gene_id = names(clust$membership), cluster = factor(clust$membership))
    betweenness <- betweenness(subnet) 
    centrality <- degree(subnet) 
    df$betweenness <- betweenness[as.character(df$gene_id)]
    df$degree_centrality <- centrality[as.character(df$gene_id)]
    df$betweenness <- as.integer(df$betweenness)
    df$degree_centrality <- as.integer(df$degree_centrality)
    df$eigen_centrality <- eigen_centrality(subnet)$vector
    df$page_rank <- page_rank(subnet)$vector
    df$closeness <- closeness(subnet)
    df$prize <- V(subnet)$prize
    df$type <- V(subnet)$type
    
    df <- df[order(-df$degree_centrality), ]
    rownames(df) <- NULL
  }
  
  cat("Done.")
  
  if (run_PCSF == FALSE) {
    return(list(DE_data = DE_data,
                PCA_data = PCA_data,
                subtype_counts_filt = counts_filt,
                sample_info = sample_info))
  } else {
    return(list(DE_data = DE_data,
                PCA_data = PCA_data,
                subtype_counts_filt = counts_filt,
                sample_info = sample_info,
                PCSF_df = df,
                PCSF_subnet = subnet))
  }
}


# Load data
load("RData/LumA/DE_data.RData")
load("RData/LumB/DE_data.RData")
load("RData/Her2/DE_data.RData")
load("RData/basal/DE_data.RData")

GTEx_data <- read.table("../BRCA_pipe/gene_reads_2017-06-05_v8_breast_mammary_tissue.gct", skip = 2)
colnames(GTEx_data) <- GTEx_data[1,]
GTEx_data <- GTEx_data[-1,-1]
rownames(GTEx_data) <- NULL

# opt for having gene Ensembl IDs instead of gene names as rownames (same as TCGA)
GTEx_ENS <- tibble::column_to_rownames(GTEx_data, "Name")
rownames(GTEx_ENS) <- gsub("\\.\\d+", "", rownames(GTEx_ENS))
GTEx_ENS <- GTEx_ENS[ , -1]
rownames <- rownames(GTEx_ENS)
GTEx_ENS <- as.data.frame(sapply(GTEx_ENS, as.numeric))
rownames(GTEx_ENS) <- rownames
rm(rownames, GTEx_data)
GTEx_ENS[] <- lapply(GTEx_ENS, function(x){as.integer(x)})
gc()

lumA_results <- subtype_analysis(subtype_data = LumA_unstranded,
                                 control_data = GTEx_ENS,
                                 subtype_name = "LumA",
                                 PCA = TRUE,
                                 run_PCSF = TRUE)

lumB_results <- subtype_analysis(subtype_data = LumB_unstranded,
                                 control_data = GTEx_ENS,
                                 subtype_name = "LumB",
                                 PCA = TRUE,
                                 run_PCSF = TRUE)

Her2_results <- subtype_analysis(subtype_data = Her2_unstranded,
                                 control_data = GTEx_ENS,
                                 subtype_name = "Her2",
                                 PCA = TRUE,
                                 run_PCSF = TRUE)

basal_results <- subtype_analysis(subtype_data = Basal_unstranded,
                                  control_data = GTEx_ENS,
                                  subtype_name = "basal",
                                  PCA = TRUE,
                                  run_PCSF = TRUE)


save(lumA_results, file = "~/OneDrive - RMIT University/PhD/large_git_files/BRCA_subtype_analysis/LumA_results.RData")
save(lumB_results, file = "~/OneDrive - RMIT University/PhD/large_git_files/BRCA_subtype_analysis/lumB_results.RData")
save(Her2_results, file = "~/OneDrive - RMIT University/PhD/large_git_files/BRCA_subtype_analysis/Her2_results.RData")
save(basal_results, file = "~/OneDrive - RMIT University/PhD/large_git_files/BRCA_subtype_analysis/basal_results.RData")

# add absolute values to df
lumA_difExp <- lumA_results$DE_data$dif_exp
lumA_difExp$logFC_abs <- abs(lumA_difExp$logFC)

lumB_difExp <- lumB_results$DE_data$dif_exp
lumB_difExp$logFC_abs <- abs(lumB_difExp$logFC)

Her2_difExp <- Her2_results$DE_data$dif_exp
Her2_difExp$logFC_abs <- abs(Her2_difExp$logFC)

basal_difExp <- basal_results$DE_data$dif_exp
basal_difExp$logFC_abs <- abs(basal_difExp$logFC)

# also save individual dif_exp results
write.csv(lumA_difExp, "latest_run/subtype_analysis/lumA_difExp.csv", row.names = FALSE)
write.csv(lumB_difExp, "latest_run/subtype_analysis/lumB_difExp.csv", row.names = FALSE)
write.csv(Her2_difExp, "latest_run/subtype_analysis/Her2_difExp.csv", row.names = FALSE)
write.csv(basal_difExp, "latest_run/subtype_analysis/basal_difExp.csv", row.names = FALSE)





