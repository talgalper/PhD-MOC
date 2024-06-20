library(biomaRt)
library(tidyverse)
library(WGCNA)
library(edgeR)
library(DESeq2)
library(matrixStats)
library(gridExtra)
library(doParallel)
library(reshape2)
library(igraph)

nCores = 8
registerDoParallel(cores = nCores)
enableWGCNAThreads(nThreads = nCores)
WGCNAnThreads()

# perform quality control
filter_low_expr <- function(tumour_matrix, control_matrix) {
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


# function to pick soft thresholding power
pick_power <- function(WGCNA_data, network_type) {
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

# load in data
load("../BRCA_pipe/RData/TCGA_normal.RData")
load("../BRCA_pipe/RData/LumA/DE_data.RData")
load("../BRCA_pipe/RData/LumB/DE_data.RData")
load("../BRCA_pipe/RData/Her2/DE_data.RData")
load("../BRCA_pipe/RData/basal/DE_data.RData")

# combine all tumour samples
all_subtypes <- cbind(LumA_unstranded, LumB_unstranded, Her2_unstranded, Basal_unstranded)

# QC + combines tumour and control samples
all_subtype_counts_filt <- filter_low_expr(tumour_matrix = all_subtypes,
                                           control_matrix = normal_unstranded)

# normalisation
all_wgcna_data <- as.matrix(all_subtype_counts_filt)
all_wgcna_data <- varianceStabilizingTransformation(all_wgcna_data)
all_wgcna_data <- as.data.frame(all_wgcna_data)
all_wgcna_data <- t(all_wgcna_data)

# clean env
rm(LumA_unstranded, LumB_unstranded, Her2_unstranded, Basal_unstranded, normal_unstranded, all_subtype_counts_filt)
collectGarbage()

# choose soft thresholding power
unsigned_sft_data <- pick_power(WGCNA_data = all_wgcna_data,
                                network_type = "unsigned")
signed_sft_data <- pick_power(WGCNA_data = all_wgcna_data,
                              network_type = "signed")

# RESTART R AND LOAD WGCNA ONLY
# identify modules: TOMType = "signed", networkType = "unsigned"
all_subtype_bwnet <- network_modules(WGCNA_data = all_wgcna_data,
                                     Power = 10)

save(all_subtype_bwnet, file = "BRCA/RData/all_TCGA/all_subtype_bwnet.RData")




# create tumour and control adj matrix
all_adjacencies <- sep_adj_matrix(WGCNA_data = all_wgcna_data,
                                  tumour_expr_df = all_subtypes,
                                  control_expr_df = normal_unstranded,
                                  power = 10)

# save adj matrix
save(all_adjacencies, file = "../../../../Desktop/WGCNA_BRCA_large_files/all_subtype_adj.RData")
load("../../../../Desktop/WGCNA_BRCA_large_files/all_subtype_adj.RData")












multidata <- multiData(Reference = benign_adj, 
                       Test = disease_adj)


multicolour <- list(Reference = bwnet$colors)








# function for diff_i method
# https://academic.oup.com/bioinformatics/article/36/9/2821/5711285?login=false
diff_i <- function(tumour_adj, control_adj) {
  sum_matrix <- tumour_adj + control_adj
  max_scores <- apply(sum_matrix, 2, max)
  normalised_scores <- sweep(sum_matrix, 2, max_scores, FUN = "/")  
  median <- median(normalised_scores)
  differential_weights <- normalised_scores - median
  
  return(differential_weights)
}

all_subtype_dif_net <- diff_i(tumour_adj = all_adjacencies$tumour,
                              control_adj = all_adjacencies$control)

# clear some memory 
rm(all_adjacencies)
collectGarbage()

all_subtype_edgeList <- melt(all_subtype_dif_net)
colnames(all_subtype_edgeList) <- c("node_1", "node_2", "weight")

edgeList_filt <- all_subtype_edgeList[all_subtype_edgeList$weight > -0.01 & all_subtype_edgeList$weight < 0.01, ]















