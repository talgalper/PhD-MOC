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

# combine all tumour samples
all_subtypes <- cbind(LumA_unstranded, LumB_unstranded, Her2_unstranded, Basal_unstranded)

# QC + combines tumour and control samples
counts_filt <- filter_low_expr(tumour_matrix = all_subtypes,
                               control_matrix = GTEx_ENS,
                               sep = T)

# normalisation (transposes matrix)
tumour_data <- vst_norm(counts_df = counts_filt$tumour)
control_data <- vst_norm(counts_df = counts_filt$control)

save(control_data, tumour_data, file = "BRCA/RData/GTEx/data_norm_filt.RData")

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
rm(lumA_info, lumB_info, her2_info, basal_info)

PCA_tumour <- plot_PCA(expr_data = tumour_data,
                       sample_info = tumour_info,
                       plot_tree = T,
                       output_plot_data = T)

PCA_control <- plot_PCA(expr_data = control_data,
                        sample_info = control_info,
                        plot_tree = T,
                        output_plot_data = T)

# identify why it looks like theres two clusters in GTEx data
#xena_GTEx_sample_info <- read_tsv("../../../../Downloads/denseDataOnlyDownload.tsv")
#table(xena_GTEx_sample_info$sample %in% colnames(GTEx_ENS))
#control_info <- merge(control_info, xena_GTEx_sample_info, by = "sample", all.x = T)
#control_info <- subset(control_info, select = c("sample", "_gender"))
#control_info[is.na(control_info)] <- "unknown"
#colnames(control_info)[2] <- "group"
#
#PCA_control <- plot_PCA(expr_data = control_data,
#                        sample_info = control_info,
#                        plot_tree = F,
#                        output_plot_data = T)


# identify outliers
#dynamicCut <- cutreeDynamic(PCA_tumour$htree, method = "tree", deepSplit = 2, pamRespectsDendro = FALSE)
#outlierSamples <- which(dynamicCut == 0)
#tumour_cleanExprData <- tumour_data[-outlierSamples, ]
#tumour_outlierSamples <- rownames(tumour_data)[outlierSamples]
#tumour_outlierSamples <- tumour_info[tumour_info$sample %in% tumour_outlierSamples, ]
#tumour_info_filt <- tumour_info[tumour_info$sample %in% rownames(tumour_cleanExprData), ]
#
#dynamicCut <- cutreeDynamic(PCA_control$htree, method = "tree", deepSplit = 2, pamRespectsDendro = FALSE)
#outlierSamples <- which(dynamicCut == 0)
#control_cleanExprData <- control_data[-outlierSamples, ]
#control_outlierSamples <- rownames(control_data)[outlierSamples]
#control_outlierSamples <- control_info[control_info$sample %in% control_outlierSamples, ]
#control_info_filt <- control_info[control_info$sample %in% rownames(control_cleanExprData), ]
#
#rm(dynamicCut, outlierSamples)
#
## re-plot PCA with outliers removed
#PCA_filt_tumour <- plot_PCA(expr_data = tumour_cleanExprData,
#                            sample_info = tumour_info_filt,
#                            plot_tree = T,
#                            output_plot_data = T)
#
#PCA_filt_control <- plot_PCA(expr_data = control_cleanExprData,
#                             sample_info = control_info_filt,
#                             plot_tree = T,
#                             output_plot_data = T)



# choose soft thresholding power
tumour_sft_data <- pick_power(WGCNA_data = tumour_data,
                                network_type = "unsigned")

control_sft_data <- pick_power(WGCNA_data = control_data,
                                     network_type = "unsigned")

# clean env
rm(normal_unstranded, LumA_unstranded, LumB_unstranded, Her2_unstranded, Basal_unstranded, all_subtypes,
   GTEx_ENS)
collectGarbage()

# RESTART R AND LOAD WGCNA ONLY
library(WGCNA)
library(doParallel)
nCores = 8
registerDoParallel(cores = nCores)
enableWGCNAThreads(nThreads = nCores)
WGCNAnThreads()

control_bwnet <- network_modules(WGCNA_data = control_data,
                                 Power = 6)
tumour_bwnet <- network_modules(WGCNA_data = tumour_data,
                                Power = 6)

save(control_bwnet, file = "BRCA/RData/GTEx/control_bwnet.RData")
save(tumour_bwnet, file = "BRCA/RData/GTEx/tumour_bwnet.RData")


# create tumour and control adj matrix
start_time <- Sys.time()

tumour_adj <- adjacency(tumour_data, power = 6, type = "unsigned")
control_adj <- adjacency(control_data, power = 6, type = "unsigned")

time_elapsed <- Sys.time() - start_time
print(time_elapsed)
rm(start_time)

# save adj matrix
save(tumour_adj, control_adj, file = "../../../../Desktop/WGCNA_BRCA_large_files/GTEx_tumour_sep_adj.RData")



# module preservation using expr data
common_genes <- intersect(colnames(tumour_data), colnames(control_data))
tumour_common <- tumour_data[, colnames(tumour_data) %in% common_genes]
control_common <- control_data[, colnames(control_data) %in% common_genes]

tumour_common_colours <- tumour_bwnet$colors[names(tumour_bwnet$colors) %in% common_genes]
control_common_colours <- control_bwnet$colors[names(control_bwnet$colors) %in% common_genes]

multidata <- multiData(Control = control_common, 
                       Tumour = tumour_common)
multicolour <- list(Control = control_common_colours,
                    Tumour = tumour_common_colours)

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
                                        networkType = "unsigned",
                                        quickCor = 1,
                                        randomSeed = 1234,
                                        verbose = 3,
                                        nPermutations = 50,
                                        referenceNetworks = 1,
                                        maxModuleSize = max(max(table(tumour_bwnet$colors)), 
                                                            max(table(control_bwnet$colors))),
                                        calculateClusterCoeff = F,
                                        parallelCalculation = T)
end_time <- Sys.time()
end_time - start_time

# plot results
modulePreservation_plt <- plot_preserved_modules(preserved_modules)

save(preserved_modules, modulePreservation_plt, file = "BRCA/RData/GTEx/GTEx_tumour_modulePreservation(n=50).RData")

# non-preserved modules
plot_data <- modulePreservation_plt$plot_data$plot_data
non_preserved_modules <- plot_data[plot_data$medianRank.pres < 8 & plot_data$Zsummary.pres > 10, ]

nonPreservedGenes <- common_genes[tumour_common_colours %in% non_preserved_modules$cluster]

temp <- approved_openTargets[approved_openTargets$`Target ID` %in% nonPreservedGenes, ]

table(unique(approved_openTargets$`Target Approved Symbol`) %in% temp$`Target Approved Symbol`)
table(unique(OpenTargets$Target.ID) %in% temp$`Target Approved Symbol`)


# for new session, re-load data
load("BRCA/RData/GTEx/GTEx_tumour_modulePreservation(n=50).RData")
load("BRCA/RData/GTEx/tumour_bwnet.RData")
load("BRCA/RData/GTEx/control_bwnet.RData")
load("BRCA/RData/GTEx/data_norm_filt.RData")

# replot
library(gridExtra)
library(ggrepel)
grid.arrange(modulePreservation_plt$meadianRank_plt, modulePreservation_plt$Zsummary_plt, ncol = 2)

# query module data
colours <- labels2colors(tumour_bwnet$colors)
tumour_kWithin <- intramodularConnectivity.fromExpr(tumour_data, colours, power = 6)
rownames(tumour_kWithin) <- colnames(tumour_data)

save(tumour_kWithin, file = "BRCA/RData/GTEx/tumour_kWithin.RData")

# get top 10 genes for connectivity for each non-preserved module
tumour_topGenes = list()
for (module in non_preserved_modules$cluster) {
  moduleGenes = names(tumour_bwnet$colors)[tumour_bwnet$colors == module]
  moduleKWithin = tumour_kWithin[moduleGenes, ]
  topModuleGenes = head(order(moduleKWithin$kWithin, decreasing = TRUE), 10)
  tumour_topGenes[[module]] = moduleGenes[topModuleGenes]
  
  rm(moduleGenes, module, topModuleGenes)
}


tumour_topGenes <- melt(tumour_topGenes)
colnames(tumour_topGenes) <- c("ensembl_id", "module")

temp <- approved_openTargets[approved_openTargets$`Target ID` %in% tumour_topGenes$ensembl_id, ]


#library(biomaRt)
#ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
#topGenes_converted <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), 
#                         filters = "ensembl_gene_id", 
#                         values = tumour_topGenes$ensembl_id, 
#                         mart = ensembl)


OpenTargets <- read.csv("../BRCA_pipe/OpenTargets_data/OpenTargets_unique_drug.csv", row.names = 1)
OpenTargets_raw <- read_tsv("../BRCA_pipe/OpenTargets_data/breast_carcinoma_known_drugs.tsv")

NIH_targets <- read.table("../BRCA_pipe/NIH_BRCA_approved_drugs.txt", sep = "\t")
colnames(NIH_targets)[1] <- "approved_drugs"
NIH_targets$approved_drugs <- toupper(NIH_targets$approved_drugs)

approved_openTargets <- merge(NIH_targets, OpenTargets_raw, by.x = "approved_drugs", by.y = "Drug Name")
approved_openTargets <- approved_openTargets[!duplicated(approved_openTargets$approved_drugs) | !duplicated(approved_openTargets$`Target ID`), ]


table(unique(approved_openTargets$`Target ID`) %in% tumour_topGenes$ensembl_id)
table(unique(OpenTargets$Target.ID) %in% tumour_topGenes$ensembl_id)

# FDA targets
temp <- tumour_kWithin[rownames(tumour_kWithin) %in% approved_openTargets$`Target ID`, ]
temp <- merge(tumour_kWithin, approved_openTargets, by.x = "row.names", by.y = "Target ID")
temp <- temp[!duplicated(temp$Row.names), ]
temp_non <- approved_openTargets[!approved_openTargets$`Target ID` %in% rownames(tumour_kWithin), ]

# OpenTargets 
temp <- tumour_kWithin[rownames(tumour_kWithin) %in% OpenTargets$Target.ID, ]
temp <- merge(tumour_kWithin, OpenTargets, by.x = "row.names", by.y = "Target.ID")
temp <- temp[!duplicated(temp$Row.names), ]


## create network
# load in data for Ubuntu
load("../../../../Desktop/WGCNA_BRCA_large_files/tumour_TOM.RData")
load("../../../../OneDrive - RMIT University/PhD/large_git_files/WGCNA/tumour_TOM.RData")
load("../../../../OneDrive - RMIT University/PhD/large_git_files/WGCNA/data_norm_filt_GTEx.RData")

# calculate TOM similarity
tumour_TOM <- TOMsimilarityFromExpr(tumour_data, power = 6, nThreads = 8)
backup_tmuour_TOM <- tumour_TOM

geneNames <- colnames(tumour_data)


# subset TOM for non-preserved genes
dimnames(tumour_TOM) = list(geneNames, geneNames)
nonPreservedGenes <- geneNames[tumour_bwnet$colors %in% non_preserved_modules$cluster]
tumour_TOM_subset <- tumour_TOM[nonPreservedGenes, nonPreservedGenes]

# Converting TOM matrix to edge list
tumour_TOM_subset[upper.tri(tumour_TOM_subset, diag = TRUE)] = NA
edgeList <- which(!is.na(tumour_TOM_subset), arr.ind = TRUE)
weights <- tumour_TOM_subset[!is.na(tumour_TOM_subset)]

# Create a data frame for igraph
edges <-  data.frame(
  from = geneNames[edgeList[, 1]],
  to = geneNames[edgeList[, 2]],
  weight = weights
)

# load graph data for Ubuntu
load("../../../../Desktop/WGCNA_BRCA_large_files/igraph_data.RData")
# load data fro mac
load("../../../../OneDrive - RMIT University/PhD/large_git_files/WGCNA/igraph_data.RData")


# subset genes from non-preserved modules



library(igraph)
# Create an igraph object
network = graph_from_data_frame(edges, directed = FALSE)

# Assign colors to modules
moduleColorsSubset = tumour_bwnet$colors[geneNames %in% nonPreservedGenes]
colorPalette = unique(moduleColorsSubset)
names(colorPalette) = unique(moduleColorsSubset)
vertexColors = colorPalette[moduleColorsSubset]
V(network)$color = vertexColors

# Plot the network
plot(network, vertex.size = 5, vertex.label = NA, edge.width = E(network)$weight, 
     vertex.color = V(network)$color, edge.color = E(network)$color)





#tumour_top_hubs <- chooseTopHubInEachModule(tumour_data, colorh = tumour_bwnet$colors)

tumour_module_pvals <- cor(tumour_bwnet$MEs, tumour_data, use = "p")
tumour_module_pvals <- corPvalueStudent(tumour_module_pvals, nSamples = nrow(tumour_data))
tumour_module_pvals <- t(tumour_module_pvals)
tumour_module_pvals <- as.data.frame(tumour_module_pvals)

top_pval_genes <- list()
for (module in colnames(tumour_module_pvals)) {
  sorted_genes <- tumour_module_pvals[order(tumour_module_pvals[[module]], decreasing = FALSE), ]
  top_genes <- rownames(sorted_genes)[1:10]
  top_pval_genes[[module]] <- top_genes
}

top_pval_genes <- melt(top_pval_genes)

temp <- approved_openTargets[approved_openTargets$`Target ID` %in% top_pval_genes$value, ]
table(unique(approved_openTargets$`Target Approved Symbol`) %in% unique(temp$`Target Approved Symbol`))







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















