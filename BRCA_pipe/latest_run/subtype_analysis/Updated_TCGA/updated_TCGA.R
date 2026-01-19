library(TCGAbiolinks)
library(SummarizedExperiment)
library(edgeR)

query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification"
)

query_output <- getResults(query)


subtype_query <- TCGAquery_subtype("brca")

BRCA_sample_data <- merge(query_output[, c("file_name","cases", "cases.submitter_id", "sample_type")], 
                          subtype_query[, c("patient", "pathologic_stage", "BRCA_Subtype_PAM50")],
                          by.x = "cases.submitter_id", by.y = "patient", all.x = TRUE)


GDCdownload(query, directory = "latest_run/subtype_analysis/Updated_TCGA")

BRCA_RNAseq <- GDCprepare(query, directory = "latest_run/subtype_analysis/Updated_TCGA")

# BRCA_RNAseq_CorOutliers <- TCGAanalyze_Preprocessing(BRCA_RNAseq)

BRCA_RNAseq_unstranded <- assay(BRCA_RNAseq, "unstranded")

BRCA_normal <- BRCA_RNAseq_unstranded[, colnames(BRCA_RNAseq_unstranded) %in% 
                                        BRCA_sample_data$cases[BRCA_sample_data$sample_type == "Solid Tissue Normal"]]
rownames(BRCA_normal) <- gsub("\\.\\d+", "", rownames(BRCA_normal))

BRCA_basal <- BRCA_RNAseq_unstranded[, colnames(BRCA_RNAseq_unstranded) %in% 
                                        BRCA_sample_data$cases[BRCA_sample_data$BRCA_Subtype_PAM50 == "Basal"]]
rownames(BRCA_basal) <- gsub("\\.\\d+", "", rownames(BRCA_basal))

BRCA_lumA <- BRCA_RNAseq_unstranded[, colnames(BRCA_RNAseq_unstranded) %in% 
                                        BRCA_sample_data$cases[BRCA_sample_data$BRCA_Subtype_PAM50 == "LumA"]]
rownames(BRCA_lumA) <- gsub("\\.\\d+", "", rownames(BRCA_lumA))

BRCA_lumB <- BRCA_RNAseq_unstranded[, colnames(BRCA_RNAseq_unstranded) %in% 
                                        BRCA_sample_data$cases[BRCA_sample_data$BRCA_Subtype_PAM50 == "LumB"]]
rownames(BRCA_lumB) <- gsub("\\.\\d+", "", rownames(BRCA_lumB))

BRCA_Her2 <- BRCA_RNAseq_unstranded[, colnames(BRCA_RNAseq_unstranded) %in% 
                                        BRCA_sample_data$cases[BRCA_sample_data$BRCA_Subtype_PAM50 == "Her2"]]
rownames(BRCA_Her2) <- gsub("\\.\\d+", "", rownames(BRCA_Her2))

write.csv(BRCA_normal, "~/OneDrive - RMIT University/PhD/large_git_files/BRCA_subtype_analysis/updated_TCGA/TCGA_normal_RNAseq.csv")
write.csv(BRCA_basal, "~/OneDrive - RMIT University/PhD/large_git_files/BRCA_subtype_analysis/updated_TCGA/TCGA_basal_RNAseq.csv")
write.csv(BRCA_lumA, "~/OneDrive - RMIT University/PhD/large_git_files/BRCA_subtype_analysis/updated_TCGA/TCGA_lumA_RNAseq.csv")
write.csv(BRCA_lumB, "~/OneDrive - RMIT University/PhD/large_git_files/BRCA_subtype_analysis/updated_TCGA/TCGA_lumB_RNAseq.csv")
write.csv(BRCA_Her2, "~/OneDrive - RMIT University/PhD/large_git_files/BRCA_subtype_analysis/updated_TCGA/TCGA_Her2_RNAseq.csv")

save(query, query_output, BRCA_sample_data, file = "latest_run/subtype_analysis/Updated_TCGA/query_data.RData")





## function to repeat differential expression analysis
subtype_analysis <- function(subtype_data, control_data, subtype_name, PCA = TRUE, run_PCSF = FALSE, confirm = TRUE) {
  suppressMessages(library(edgeR))
  # laod in good DE functions from MOC_pipe
  source("../MOC_pipe/R/functions.R")
  
  # create sample info
  sample_info <- data.frame(sample = c(colnames(subtype_data), colnames(control_data)),
                            group = c(rep(subtype_name, ncol(subtype_data)), rep("Control", ncol(control_data))))
  
  # combine sample type data
  load("latest_run/subtype_analysis/Updated_TCGA/query_data.RData")
  BRCA_sample_data$cases <- gsub("-", "\\.", BRCA_sample_data$cases)
  sample_info <- merge(sample_info, BRCA_sample_data[,c(3:6)], by.x = "sample", by.y = "cases", all.x = T)
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
  
  if (isTRUE(confirm)) {
    # offer abort if PCA doest look good
    ans <- readline(prompt = "Does everything look ok? continue... (y/n): ")
    if (tolower(ans) != "y") {
      cat("Aborted.\n")
      return(invisible(NULL))
    }
  }
  
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
  
  cat("\n", "Done.")
  
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


BRCA_lumA <- read.csv("~/OneDrive - RMIT University/PhD/large_git_files/BRCA_subtype_analysis/updated_TCGA/TCGA_lumA_RNAseq.csv", row.names = 1)
BRCA_lumB <- read.csv("~/OneDrive - RMIT University/PhD/large_git_files/BRCA_subtype_analysis/updated_TCGA/TCGA_lumB_RNAseq.csv", row.names = 1)
BRCA_Her2 <- read.csv("~/OneDrive - RMIT University/PhD/large_git_files/BRCA_subtype_analysis/updated_TCGA/TCGA_Her2_RNAseq.csv", row.names = 1)
BRCA_basal <- read.csv("~/OneDrive - RMIT University/PhD/large_git_files/BRCA_subtype_analysis/updated_TCGA/TCGA_basal_RNAseq.csv", row.names = 1)
BRCA_normal <- read.csv("~/OneDrive - RMIT University/PhD/large_git_files/BRCA_subtype_analysis/updated_TCGA/TCGA_normal_RNAseq.csv", row.names = 1)

# read in and format GTEx data
GTEx_data <- read.table("../BRCA_pipe/latest_run/subtype_analysis/Updated_TCGA/GTEx_breast_mammary_v10.gct", skip = 2)
colnames(GTEx_data) <- GTEx_data[1, ]
GTEx_data <- GTEx_data[-1, ]
rownames(GTEx_data) <- NULL

# opt for having gene Ensembl IDs instead of gene names as rownames (same as TCGA)
GTEx_ENS <- textshape::column_to_rownames(GTEx_data, "Name")
rownames(GTEx_ENS) <- gsub("\\.\\d+", "", rownames(GTEx_ENS))
GTEx_ENS <- GTEx_ENS[ , -1]
rownames <- rownames(GTEx_ENS)
GTEx_ENS <- as.data.frame(sapply(GTEx_ENS, as.numeric))
rownames(GTEx_ENS) <- rownames
rm(rownames, GTEx_data)
GTEx_ENS[] <- lapply(GTEx_ENS, function(x){as.integer(x)})


lumA_results <- subtype_analysis(subtype_data = BRCA_lumA,
                                 control_data = GTEx_ENS,
                                 subtype_name = "LumA",
                                 PCA = TRUE,
                                 confirm = FALSE)

lumB_results <- subtype_analysis(subtype_data = BRCA_lumB,
                                 control_data = GTEx_ENS,
                                 subtype_name = "LumB",
                                 PCA = TRUE,
                                 confirm = FALSE)

Her2_results <- subtype_analysis(subtype_data = BRCA_Her2,
                                 control_data = GTEx_ENS,
                                 subtype_name = "Her2",
                                 PCA = TRUE,
                                 confirm = FALSE)

basal_results <- subtype_analysis(subtype_data = BRCA_basal,
                                  control_data = GTEx_ENS,
                                  subtype_name = "Basal",
                                  PCA = TRUE,
                                  confirm = FALSE)



save(lumA_results, file = "~/OneDrive - RMIT University/PhD/large_git_files/BRCA_subtype_analysis/updated_TCGA/LumA_results.RData")
save(lumB_results, file = "~/OneDrive - RMIT University/PhD/large_git_files/BRCA_subtype_analysis/updated_TCGA/lumB_results.RData")
save(Her2_results, file = "~/OneDrive - RMIT University/PhD/large_git_files/BRCA_subtype_analysis/updated_TCGA/Her2_results.RData")
save(basal_results, file = "~/OneDrive - RMIT University/PhD/large_git_files/BRCA_subtype_analysis/updated_TCGA/basal_results.RData")



basal_dif_exp <- basal_results$DE_data$dif_exp
basal_dif_exp <- id_annot(ensembl, data = basal_dif_exp, 
                          input_type = "ensembl_gene_id", 
                          convert_to = "external_gene_name")

basal_difExp <- id_annot(ensembl = ensembl, 
                         data = basal_difExp,
                         input_type = "ensembl_gene_id", 
                         convert_to = "external_gene_name")





