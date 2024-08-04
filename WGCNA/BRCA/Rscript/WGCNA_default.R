### This script combines all the methods used in the different WGCNA attempts to hopefully be the best one ###
## firstly, combines TCGA and GTEx data to perform trait correlation ###
### secondly, analyses the disease groups seperately for module preservation ###

library(tidyverse)
library(WGCNA)
library(edgeR)
library(DESeq2)
library(doParallel)

nCores = 8
registerDoParallel(cores = nCores)
enableWGCNAThreads(nThreads = nCores)
WGCNAnThreads()


# load in data
load("../BRCA_pipe/RData/LumA/DE_data.RData")
load("../BRCA_pipe/RData/LumB/DE_data.RData")
load("../BRCA_pipe/RData/Her2/DE_data.RData")
load("../BRCA_pipe/RData/basal/DE_data.RData")
load("../BRCA_pipe/RData/TCGA_normal.RData")

# load data
GTEx_data <- read.table("../BRCA_pipe/gene_reads_2017-06-05_v8_breast_mammary_tissue.gct", skip = 2)
colnames(GTEx_data) <- GTEx_data[1, ]
GTEx_data <- GTEx_data[-1, -1]
rownames(GTEx_data) <- NULL

# opt for having gene Ensembl IDs instead of gene names as rownames (same as TCGA)
GTEx_ENS <- column_to_rownames(GTEx_data, "Name")
rownames(GTEx_ENS) <- gsub("\\.\\d+", "", rownames(GTEx_ENS))
GTEx_ENS <- GTEx_ENS[ , -1]
rownames <- rownames(GTEx_ENS)
GTEx_ENS <- as.data.frame(sapply(GTEx_ENS, as.numeric))
rownames(GTEx_ENS) <- rownames
rm(rownames, GTEx_data)
GTEx_ENS[] <- lapply(GTEx_ENS, function(x){as.integer(x)})

# sample info
control_info <- data.frame(sample = colnames(GTEx_ENS),
                           group = rep("control", ncol(GTEx_ENS)))
lumA_info <- data.frame(sample = colnames(LumA_unstranded),
                        group = rep("lumA", ncol(LumA_unstranded)))
lumB_info <- data.frame(sample = colnames(LumB_unstranded),
                        group = rep("lumB", ncol(LumB_unstranded)))
her2_info <- data.frame(sample = colnames(Her2_unstranded),
                        group = rep("Her2", ncol(Her2_unstranded)))
basal_info <- data.frame(sample = colnames(Basal_unstranded),
                         group = rep("basal", ncol(Basal_unstranded)))
sample_info <- rbind(control_info, lumA_info, lumB_info, her2_info, basal_info)

# combine all tumour samples
all_subtypes <- cbind(LumA_unstranded, LumB_unstranded, Her2_unstranded, Basal_unstranded)

# clean env
rm(normal_unstranded, LumA_unstranded, LumB_unstranded, Her2_unstranded, Basal_unstranded)
rm(control_info, lumA_info, lumB_info, her2_info, basal_info)
collectGarbage()

# QC + combines tumour and control samples
all_subtype_counts_filt <- filter_low_expr(tumour_matrix = all_subtypes,
                                           control_matrix = GTEx_ENS)

# normalisation (transposes matrix)
all_wgcna_data <- vst_norm(all_subtype_counts_filt)

# plot PCA
PCA_results <- plot_PCA(expr_data = all_wgcna_data,
                        sample_info = sample_info,
                        plot_tree = T,
                        output_plot_data = T)


sft <- pickSoftThreshold(all_wgcna_data,
                         blockSize = 45000,
                         verbose = 2)
sft <- sft$fitIndices

library(gridExtra)
library(grid)
a1 <- ggplot(sft, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit') +
  theme_classic()

a2 <- ggplot(sft, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()

grid.arrange(a1, a2, nrow = 2)
rm(a1, a2)

save(sft, file = "BRCA/RData/all_default/combined_sft.RData")


# RESTART R AND LOAD WGCNA ONLY
library(WGCNA)
library(doParallel)
nCores = 8
registerDoParallel(cores = nCores)
enableWGCNAThreads(nThreads = nCores)
WGCNAnThreads()

start_time <- Sys.time()
bwnet <- blockwiseModules(all_wgcna_data,
                          maxBlockSize = 45000,
                          power = 6,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3,
                          saveTOMs = FALSE)
elapsed_time <- Sys.time() - start_time
print(paste0("Elapsed time: ", elapsed_time))


# Plot the dendrogram
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)


save(bwnet, file = "BRCA/RData/all_default/all_bwnet.RData")
load("BRCA/RData/all_default/all_bwnet.RData")


## trait correlation
library(tidyverse)
library(CorLevelPlot)

sample_info <- data.frame(row.names = sample_info$sample,
                          status = c(rep("tumour", ncol(all_subtypes)), rep("control", ncol(GTEx_ENS))),
                          group = c(rep("lumA", ncol(LumA_unstranded)),
                                    rep("GTEx", ncol(GTEx_ENS)),
                                    rep("lumB", ncol(LumB_unstranded)),
                                    rep("Her2", ncol(Her2_unstranded)),
                                    rep("basal", ncol(Basal_unstranded))))



traits.state <- binarizeCategoricalColumns.forPlots(sample_info$status)
traits.subtype <- binarizeCategoricalColumns.forPlots(sample_info$group)
traits <- cbind(traits.state, traits.subtype)
rownames(traits) <- c(colnames(all_subtypes), colnames(GTEx_ENS))

moduleTrait_cor <- cor(bwnet$MEs, traits, use = "p")
moduleTrait_cor_pvals <- corPvalueStudent(moduleTrait_cor, nSamples = nrow(all_wgcna_data))

heatmap_data <- merge(bwnet$MEs, traits, by = "row.names")

heatmap_data <- column_to_rownames(heatmap_data, "Row.names")

CorLevelPlot(heatmap_data,
             x = names(heatmap_data)[18:24],
             y = names(heatmap_data)[1:17],
             col = c("blue1", "skyblue", "white", "pink", "red"))




# perform GO and pathway analysis





###############################################################################
# Will now perform WGCNA on TCGA and GTEx groups separately for module     
# preservation analysis                                                    
###############################################################################
library(tidyverse)
library(WGCNA)
library(edgeR)
library(DESeq2)
library(matrixStats)
library(doParallel)
library(reshape2)
library(igraph)

nCores = 8
registerDoParallel(cores = nCores)
enableWGCNAThreads(nThreads = nCores)
WGCNAnThreads()

# load in data
load("../BRCA_pipe/RData/LumA/DE_data.RData")
load("../BRCA_pipe/RData/LumB/DE_data.RData")
load("../BRCA_pipe/RData/Her2/DE_data.RData")
load("../BRCA_pipe/RData/basal/DE_data.RData")
load("../BRCA_pipe/RData/TCGA_normal.RData")

# load data
GTEx_data <- read.table("../BRCA_pipe/gene_reads_2017-06-05_v8_breast_mammary_tissue.gct", skip = 2)
colnames(GTEx_data) <- GTEx_data[1, ]
GTEx_data <- GTEx_data[-1, -1]
rownames(GTEx_data) <- NULL

# opt for having gene Ensembl IDs instead of gene names as rownames (same as TCGA)
GTEx_ENS <- column_to_rownames(GTEx_data, "Name")
rownames(GTEx_ENS) <- gsub("\\.\\d+", "", rownames(GTEx_ENS))
GTEx_ENS <- GTEx_ENS[ , -1]
rownames <- rownames(GTEx_ENS)
GTEx_ENS <- as.data.frame(sapply(GTEx_ENS, as.numeric))
rownames(GTEx_ENS) <- rownames
rm(rownames, GTEx_data)
GTEx_ENS[] <- lapply(GTEx_ENS, function(x){as.integer(x)})

# load in QC data
load("../../../../Desktop/WGCNA_BRCA_large_files/data_norm_filt_GTEx.RData")

# plot PCA
control_info <- data.frame(sample = colnames(GTEx_ENS),
                           group = rep("control", ncol(GTEx_ENS)))

lumA_info <- data.frame(sample = colnames(LumA_unstranded),
                        group = rep("lumA", ncol(LumA_unstranded)))
lumB_info <- data.frame(sample = colnames(LumB_unstranded),
                        group = rep("lumB", ncol(LumB_unstranded)))
her2_info <- data.frame(sample = colnames(Her2_unstranded),
                        group = rep("Her2", ncol(Her2_unstranded)))
basal_info <- data.frame(sample = colnames(Basal_unstranded),
                         group = rep("basal", ncol(Basal_unstranded)))
tumour_info <- rbind(lumA_info, lumB_info, her2_info, basal_info)

# clean env
rm(normal_unstranded, LumA_unstranded, LumB_unstranded, Her2_unstranded, Basal_unstranded)
rm(lumA_info, lumB_info, her2_info, basal_info)
collectGarbage()

# plot PCA
PCA_tumour <- plot_PCA(expr_data = tumour_data,
                       sample_info = tumour_info,
                       plot_tree = T,
                       output_plot_data = T)

PCA_control <- plot_PCA(expr_data = control_data,
                        sample_info = control_info,
                        plot_tree = T,
                        output_plot_data = T)


# choose soft thresholding power
pick_power <- function(WGCNA_data) {
  sft <- pickSoftThreshold(WGCNA_data,
                           blockSize = 45000,
                           verbose = 2)
  sft <- sft$fitIndices
  
  library(gridExtra)
  library(grid)
  a1 <- ggplot(sft, aes(Power, SFT.R.sq, label = Power)) +
    geom_point() +
    geom_text(nudge_y = 0.1) +
    geom_hline(yintercept = 0.8, color = 'red') +
    labs(x = 'Power', y = 'Scale free topology model fit') +
    theme_classic()
  
  a2 <- ggplot(sft, aes(Power, mean.k., label = Power)) +
    geom_point() +
    geom_text(nudge_y = 0.1) +
    labs(x = 'Power', y = 'Mean Connectivity') +
    theme_classic()
  
  grid.arrange(a1, a2, nrow = 2, top = textGrob(deparse(substitute(WGCNA_data))))
  
  return(sft)
}

tumour_sft <- pick_power(tumour_data)
control_sft <- pick_power(control_data)

save(tumour_sft, file = "BRCA/RData/all_default/tumour_sft.RData")
save(control_sft, file = "BRCA/RData/all_default/control_sft.RData")


# RESTART R AND LOAD WGCNA ONLY
library(WGCNA)
library(doParallel)
nCores = 8
registerDoParallel(cores = nCores)
enableWGCNAThreads(nThreads = nCores)
WGCNAnThreads()

network_modules <- function(WGCNA_data, Power) {
  start_time <- Sys.time()
  bwnet <- blockwiseModules(WGCNA_data,
                            maxBlockSize = 45000,
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

tumour_bwnet <- network_modules(tumour_data, Power = 6)
control_bwnet <- network_modules(control_data, Power = 6)

save(tumour_bwnet, file = "BRCA/RData/all_default/tumour_bwnet.RData")
save(control_bwnet, file = "BRCA/RData/all_default/control_bwnet.RData")

load("BRCA/RData/all_default/tumour_bwnet.RData")
load("BRCA/RData/all_default/control_bwnet.RData")


# module preservation analysis
multidata <- multiData(Control = control_data, 
                       Tumour = tumour_data)
multicolour <- list(Control = control_bwnet$colors,
                    Tumour = tumour_bwnet$colors)

# RESTART R AND LOAD WGCNA ONLY
library(WGCNA)
library(doParallel)
nCores = 8
registerDoParallel(cores = nCores)
enableWGCNAThreads(nThreads = nCores)
WGCNAnThreads()

start_time <- Sys.time()
preserved_modules <- modulePreservation(multiData = multidata,
                                        multiColor = multicolour,
                                        dataIsExpr = T,
                                        quickCor = 1,
                                        randomSeed = 1234,
                                        verbose = 3,
                                        nPermutations = 100,
                                        maxModuleSize = max(max(table(tumour_bwnet$colors)), 
                                                            max(table(control_bwnet$colors))),                                        calculateClusterCoeff = F,
                                        parallelCalculation = T)
end_time <- Sys.time()
end_time - start_time

# plot results
modulePreservation_plt <- plot_preserved_modules(preserved_modules)

save(preserved_modules, modulePreservation_plt, file = "BRCA/RData/all_default/modulePreservation(n=100).RData")

# non preserved modules
plot_data <- modulePreservation_plt$plot_data$plot_data
non_preserved_modules <- plot_data[plot_data$medianRank.pres > 8 & plot_data$Zsummary.pres < 10, ]

# plot cross-tabulation
cross_tab_counts <- preserved_modules$accuracy$observedCounts$ref.Control$inColumnsAlsoPresentIn.Tumour
cross_tab_pvalues <- preserved_modules$accuracy$observedFisherPvalues$ref.Control$inColumnsAlsoPresentIn.Tumour

data <- melt(cross_tab_counts)
colnames(data) <- c("Test", "Reference", "Count")

pvalues_melted <- melt(cross_tab_pvalues)
colnames(pvalues_melted) <- c("Test", "Reference", "PValue")
# Combine the counts and p-values
data$PValue <- pvalues_melted$PValue
data$LogPValue <- -log10(data$PValue)

# Add the total size of each module
control_sizes <- as.data.frame(table(control_bwnet$colors))
colnames(control_sizes) <- c("Module", "ControlSize")

tumour_sizes <- as.data.frame(table(tumour_bwnet$colors))
colnames(tumour_sizes) <- c("Module", "TumourSize")

# Merge the sizes with the data
data <- merge(data, tumour_sizes, by.x = "Test", by.y = "Module")
data <- merge(data, control_sizes, by.x = "Reference", by.y = "Module")

# Create labels with module size
data$Reference <- paste(data$Reference, "(", data$ControlSize, ")", sep = "")
data$Test <- paste(data$Test, "(", data$TumourSize, ")", sep = "")

ggplot(data, aes(x = Test, y = Reference, fill = LogPValue)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "red") +
  geom_text(aes(label = Count), size = 3) +  # Display only counts
  theme_minimal() +
  labs(title = "Control modules (rows) vs. Tumour modules (columns)",
       x = "Test Modules", y = "Reference Modules") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# perform enrichment on modules
library(clusterProfiler)
library(org.Hs.eg.db)
library(progress)
pb <- progress_bar$new(total = length(unique(control_bwnet$colors)))

# run GO enrichment on control modules
control_GO <- list()
for (module in unique(control_bwnet$colors)) {
  genes <- names(control_bwnet$colors)[control_bwnet$colors %in% module]
  GO <- enrichGO(genes, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "ALL")
  result <- GO@result
  split_result <- split(result, result$ONTOLOGY)
  result_first_5 <- lapply(split_result, function(df) head(df, 5))
  result_first_5 <- do.call(rbind, result_first_5)
  
  control_GO[[module]] <- result_first_5
  
  rm(genes, GO, split_result, result_first_5, module)
  pb$tick()
}

# run GO enrichment on tumour modules
pb <- progress_bar$new(total = length(unique(tumour_bwnet$colors)))
tumour_GO <- list()
for (module in unique(tumour_bwnet$colors)) {
  genes <- names(tumour_bwnet$colors)[tumour_bwnet$colors %in% module]
  GO <- enrichGO(genes, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "ALL")
  result <- GO@result
  split_result <- split(result, result$ONTOLOGY)
  result_first_5 <- lapply(split_result, function(df) head(df, 5))
  result_first_5 <- do.call(rbind, result_first_5)
  
  tumour_GO[[module]] <- result_first_5
  
  rm(genes, GO, split_result, result_first_5, module)
  pb$tick()
}
rm(pb)

save(control_GO, tumour_GO, file = "BRCA/RData/all_default/module_GO_data.RData")








