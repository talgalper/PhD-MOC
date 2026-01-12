### Updated analysis of individual BRCA subtypes ###
### This analysis is in response to thesis reveiwer feedback ###

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



### Rank sensitivity analysis ###

# read in data
library(progress)
library(biomaRt)
library(data.table)
source("../MOC_pipe/R/functions.R")

druggability <- read.csv("../Druggability_analysis/data_general/druggability_scores_annot.csv")
load("../ML/RData/full_fpocket_results.RData")
results_master$Score <- (results_master$Score - min(results_master$Score)) / (max(results_master$Score) - min(results_master$Score))
results_master <- results_master[order(-results_master$`Druggability Score`, -results_master$Score), ]
results_master <- results_master[!duplicated(results_master$uniprot_id), ]
druggability <- subset(druggability, 
                       select = c("external_gene_name", "druggability", "CP_score", 
                                  "highest_score", "source"))
druggability <- unique(druggability)

# PubTator3 counts and get total counts per gene
PubTator <- fread("~/OneDrive - RMIT University/PhD/large_git_files/PubTator3/citaiton_counts_ognsmAnnot.csv")
PubTator$tax_id <- as.character(PubTator$tax_id)
PubTator$entrezgene_id <- as.character(PubTator$entrezgene_id)
PubTator[, combined := ifelse(symbol == "", entrezgene_id, symbol)]
PubTator_counts <- PubTator[
  , .(counts = sum(count)),
  by = combined
]
PubTator_counts <- PubTator_counts[order(-PubTator_counts$counts), ]
colnames(PubTator_counts)[1] <- "symbol"

# load in lightweight ML data
ML_scores <- read.csv("../ML/results/ML_scores.csv")

## Rank sensitivity function
avg_rank_sensitivity <- function(input_data,
                                 features = list(
                                   betweenness = "betweenness",
                                   centrality = "degree_centrality",
                                   druggability = "highest_score",
                                   eigen_centrality = "eigen_centrality",
                                   closeness = "closeness",
                                   citation = "counts_norm"
                                 ),
                                 # Specify which features should be subtracted (e.g. citation counts)
                                 negative_features = c("citation"),
                                 step = 0.1) {
  # Extract feature names and build corresponding weight column names.
  feature_names <- names(features)
  weight_names <- paste0(feature_names, "_w")
  
  # Dynamically generate all weight combinations for the selected features.
  grid_list <- lapply(seq_along(feature_names), function(i) seq(0, 1, step))
  names(grid_list) <- weight_names
  all_combos <- expand.grid(grid_list)
  
  # Keep only the combinations where the weights sum to 1 (within a tiny tolerance).
  valid_combos <- subset(all_combos, abs(rowSums(all_combos) - 1) < 1e-9)
  
  num_combos <- nrow(valid_combos)
  num_genes <- nrow(input_data)
  gene_names <- input_data$external_gene_name
  
  # Create a matrix to store the rank of every gene for each weight combination.
  rank_matrix <- matrix(NA, nrow = num_genes, ncol = num_combos)
  rownames(rank_matrix) <- gene_names
  colnames(rank_matrix) <- paste0("Combo_", seq_len(num_combos))
  
  # Set up a progress bar
  pb <- progress_bar$new(
    format = "[:bar] :current/:total (:percent) eta: :eta", 
    total = num_combos
  )
  
  # Loop over each valid weight combination.
  for(i in seq_len(num_combos)) {
    pb$tick()
    current_weights <- valid_combos[i, ]
    
    # Compute combined score dynamically based on the selected features.
    # For each feature, multiply the weight by the corresponding input_data column.
    # Use a negative sign for features specified in negative_features.
    combined_score <- rep(0, num_genes)
    for(feature in feature_names) {
      weight_col <- paste0(feature, "_w")
      sign_factor <- ifelse(feature %in% negative_features, -1, 1)
      combined_score <- combined_score +
        sign_factor * current_weights[[weight_col]] * input_data[[ features[[feature]] ]]
    }
    
    # Rank all genes: genes with higher combined scores receive lower (better) ranks.
    order_indices <- order(combined_score, decreasing = TRUE)
    ranks <- integer(num_genes)
    ranks[order_indices] <- seq_along(combined_score)
    
    # Store these ranks in the rank matrix.
    rank_matrix[, i] <- ranks
  }
  
  # Calculate the average rank and variance for each gene over all weight combinations.
  avg_ranks <- rowMeans(rank_matrix)
  var_ranks <- apply(rank_matrix, 1, var)
  
  # Create a summary data frame.
  result_df <- data.frame(
    gene = gene_names,
    avg_rank = avg_ranks,
    rank_variance = var_ranks,
    stringsAsFactors = FALSE
  )
  
  # Order the result by average rank (lower rank = higher priority).
  result_df <- result_df[order(result_df$avg_rank), ]
  
  return(result_df)
}



# read in HHnet results and loop over each BRCA subtype
load("../Hierarchical_HotNet/BRCA/subtype_analysis/subtype_results.RData")
HHnet_results <- list(lumA_results = lumA_results, 
                      lumB_results = lumB_results, 
                      Her2_results = Her2_results, 
                      basal_results = basal_results)
rm(lumA_results, lumB_results, Her2_results, basal_results)


RS_HHnet_results <- list()
for (i in seq_along(HHnet_results)) {
  # initialise subtype varoiables
  name <- names(HHnet_results[i])
  cat("Analysising for:", name, "\n")
  subtype <- HHnet_results[[i]]
  
  HHnet_df <- subtype$df_subnetNeighs
  colnames(HHnet_df)[5] <- "degree_centrality"
  
  # merge HHnet data with druggability and pubtator data
  HHnet_df <- merge.data.table(as.data.table(HHnet_df), PubTator_counts, by.x = "external_gene_name", by.y = "symbol", all.x = T)
  HHnet_df <- HHnet_df[order(-HHnet_df$degree), ]
  HHnet_df <- merge(HHnet_df, druggability, by = "external_gene_name", all.x = T)
  temp <- HHnet_df[is.na(HHnet_df$druggability), ]
  HHnet_df <- HHnet_df[complete.cases(HHnet_df), ]
  
  # number of items with missing entries
  cat("Incomplete cases:", nrow(temp), "\n")
  
  # normalise data
  HHnet_df$counts_norm <- log10(HHnet_df$counts)
  HHnet_df[, c(5:9,15)] <- lapply(HHnet_df[, c(5:9,15)], function(x) {
    (x - min(x)) / (max(x) - min(x))
  })
  
  # run RS function without citation score
  cat("Running RS algorith:", "\n")
  RS_HHnet_enrich <- avg_rank_sensitivity(HHnet_df,
                                          features = list(
                                            betweenness = "betweenness",
                                            centrality = "degree_centrality",
                                            druggability = "highest_score",
                                            eigen_centrality = "eigen_centrality",
                                            closeness = "closeness",
                                            page_rank = "page_rank"), 
                                          step = 0.1)
  
  rownames(RS_HHnet_enrich) <- NULL
  
  # annotate with uniprot ids
  cat("Annontating with uniprot IDs...")

  RS_HHnet_enrich <- id_annot(ensembl, RS_HHnet_enrich, input_type = "external_gene_name", convert_to = c("uniprot_gn_id", "description"))

  # annotate with citation counts, ML scores and druggablity
  RS_HHnet_enrich <- merge(RS_HHnet_enrich, ML_scores, by.x = "uniprot_gn_id", by.y = "Protein", all.x = T)
  RS_HHnet_enrich <- merge.data.table(RS_HHnet_enrich, PubTator_counts, by.x = "external_gene_name", by.y = "symbol", all.x = T)
  RS_HHnet_enrich <- merge(RS_HHnet_enrich, druggability, by = c("external_gene_name"), all.x = TRUE)
  
  # format cols
  RS_HHnet_enrich$CP_score <- round(RS_HHnet_enrich$CP_score, digits = 3)
  RS_HHnet_enrich$`FP/CP` <- paste(RS_HHnet_enrich$druggability, RS_HHnet_enrich$CP_score, sep = "/")
  RS_HHnet_enrich <- subset(RS_HHnet_enrich, 
                            select = c("external_gene_name", "uniprot_gn_id", "description", 
                                       "avg_rank", "rank_variance", "Prediction_Score_rf", 
                                       "counts", "FP/CP", "source"))
  
  # reorder
  RS_HHnet_enrich <- RS_HHnet_enrich[order(RS_HHnet_enrich$avg_rank), ]
  rownames(RS_HHnet_enrich) <- NULL
  
  RS_HHnet_results[[name]] <- RS_HHnet_enrich
}

rm(subtype,temp,i,name,HHnet_df,RS_HHnet_enrich)
gc()

save(RS_HHnet_results, file = "latest_run/subtype_analysis/RS_HHnet_subtype_results.RData")




load("latest_run/subtype_analysis/RS_HHnet_subtype_results.RData")


# get BRCA specific subtypes
TTD_approved_cancer <- read.csv("../Druggability_analysis/data_general/TTD_approved_EXT.csv")

# Her2 subtype
filter_terms <- c("Trastuzumab", "Pertuzumab", "Neratinib", "lapatinib")
pattern <- paste(filter_terms, collapse = "|")
TTD_approved_cancer <- TTD_approved[grepl(pattern, TTD_approved$INDICATION, ignore.case = TRUE), ]


Her2 <- RS_HHnet_results$Her2_results

# subset genes with less than 5000 citations
filt <- Her2[Her2$counts < 5000, ]
# predicted by RF
filt <- Her2[Her2$Prediction_Score_rf > 0.5, ]
rownames(filt) <- NULL

















