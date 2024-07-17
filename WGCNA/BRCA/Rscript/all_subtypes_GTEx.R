library(biomaRt)
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

rm(normal_unstranded, LumA_unstranded, LumB_unstranded, Her2_unstranded, Basal_unstranded)

# QC + combines tumour and control samples
all_subtype_counts_filt <- filter_low_expr(tumour_matrix = all_subtypes,
                                           control_matrix = GTEx_ENS)

# normalisation (transposes matrix)
all_wgcna_data <- vst_norm(all_subtype_counts_filt)

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
sample_info <- rbind(control_info, lumA_info, lumB_info, her2_info, basal_info)

PCA_results_GTEx <- plot_PCA(expr_data = all_wgcna_data,
                             sample_info = sample_info,
                             plot_tree = T,
                             output_plot_data = T)


# experimented with logCPM normalisation instead of vst
#expr_data <- cpm(as.matrix(all_subtype_counts_filt), log = T)
#PCA_results_GTEx <- plot_PCA(expr_data = t(expr_data),
#                             sample_info = sample_info,
#                             plot_tree = F,
#                             output_plot_data = T)



# identify outliers
#dynamicCut <- cutreeDynamic(PCA_results_GTEx$htree, distM = dist(all_wgcna_data), method = "tree", deepSplit = 2, pamRespectsDendro = FALSE)
#outlierSamples <- which(dynamicCut == 0)
#cleanExprData <- all_wgcna_data[-outlierSamples, ]
#outlierSamples <- all_wgcna_data[outlierSamples, ]
#outlierSamples <- sample_info[sample_info$sample %in% rownames(outlierSamples), ]
#
#sample_info_filt <- sample_info[sample_info$sample %in% rownames(cleanExprData), ]
#
## re-plot PCA with outliers removed
#PCA_results_filt_GTEx <- plot_PCA(expr_data = cleanExprData,
#                                  sample_info = sample_info_filt,
#                                  plot_tree = T,
#                                  output_plot_data = T)
#
## choose soft thresholding power
sft_data_unsigned <- pick_power(WGCNA_data = all_wgcna_data,
                                network_type = "unsigned")

sft_data_signed <- pick_power(WGCNA_data = all_wgcna_data,
                              network_type = "signed")


#sft_cleanData_unsigned <- pick_power(WGCNA_data = cleanExprData,
#                                     network_type = "unsigned")


# identify modules: TOMType = "signed", networkType = "unsigned"
# split data
#control_expr <- cleanExprData[rownames(cleanExprData) %in% colnames(normal_unstranded), ]
#tumour_expr <- cleanExprData[!rownames(cleanExprData) %in% colnames(normal_unstranded), ]

# clean env
rm(control_info, lumA_info, lumB_info, her2_info, basal_info)
collectGarbage()

# RESTART R AND LOAD WGCNA ONLY
library(WGCNA)
library(doParallel)
nCores = 8
registerDoParallel(cores = nCores)
enableWGCNAThreads(nThreads = nCores)
WGCNAnThreads()

bwnet <- network_modules(WGCNA_data = all_wgcna_data,
                         Power = 6)


save(bwnet, file = "BRCA/RData/all_together/bwnet.RData")
load("BRCA/RData/all_together/bwnet.RData")




# module preservation using expr data
# split expr data but use colours from combined
all_subtype_counts_filt <- filter_low_expr(tumour_matrix = all_subtypes,
                                           control_matrix = GTEx_ENS,
                                           sep = T)
# normalisation (transposes matrix)
tumour_wgcna_data <- vst_norm(all_subtype_counts_filt$tumour)
control_wgcna_data <- vst_norm(all_subtype_counts_filt$control)

common_genes <- intersect(colnames(tumour_wgcna_data), colnames(control_wgcna_data))
tumour_common <- tumour_wgcna_data[, colnames(tumour_wgcna_data) %in% common_genes]
control_common <- control_wgcna_data[, colnames(control_wgcna_data) %in% common_genes]

colours <- bwnet$colors
colours_common <- colours[names(colours) %in% common_genes]

multidata <- multiData(Control = control_common, 
                       Tumour = tumour_common)
multicolour <- list(Control = colours_common)

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
                                        maxModuleSize = max(table(colours_common)),
                                        calculateClusterCoeff = F,
                                        parallelCalculation = T)
end_time <- Sys.time()
end_time - start_time

# plot results
modulePreservation_plt <- plot_preserved_modules(preserved_modules)

save(preserved_modules, modulePreservation_plt, file = "BRCA/RData/all_together/modulePreservation(n=50).RData")

# non-preserved modules
plot_data <- modulePreservation_plt$plot_data$plot_data
non_preserved_modules <- plot_data[plot_data$medianRank.pres > 8 & plot_data$Zsummary.pres < 10, ]

nonPreservedGenes <- common_genes[tumour_common_colours %in% non_preserved_modules$cluster]

temp <- approved_openTargets[approved_openTargets$`Target ID` %in% nonPreservedGenes, ]

table(unique(approved_openTargets$`Target Approved Symbol`) %in% temp$`Target Approved Symbol`)
table(unique(OpenTargets$Target.ID) %in% temp$`Target Approved Symbol`)



























### trait correlation
sample_info <- data.frame(row.names = c(colnames(all_subtypes), colnames(GTEx_ENS)),
                          status = c(rep("tumour", ncol(all_subtypes)), rep("control", ncol(GTEx_ENS))),
                          group = c(rep("lumA", ncol(LumA_unstranded)),
                                    rep("GTEx", ncol(GTEx_ENS)),
                                    rep("lumB", ncol(LumB_unstranded)),
                                    rep("Her2", ncol(Her2_unstranded)),
                                    rep("basal", ncol(Basal_unstranded))))


traits.state <- binarizeCategoricalColumns(sample_info$status,
                                     includePairwise = F,
                                     includeLevelVsAll = T)

traits.subtype <- binarizeCategoricalColumns(sample_info$group,
                                             includePairwise = F,
                                             includeLevelVsAll = T,
                                             levelOrder = c("GTEx", "lumA", "lumB", "Her2", "basal"))

traits <- cbind(traits.state, traits.subtype)
rownames(traits) <- c(colnames(all_subtypes), colnames(GTEx_ENS))

moduleTrait_cor <- cor(all_subtype_bwnet$MEs, traits, use = "p")
moduleTrait_cor_pvals <- corPvalueStudent(moduleTrait_cor, nSamples = nrow(all_wgcna_data))

heatmap_data <- merge(all_subtype_bwnet$MEs, traits, by = "row.names")

heatmap_data <- column_to_rownames(heatmap_data, "Row.names")

CorLevelPlot(heatmap_data,
             x = names(heatmap_data)[17:21],
             y = names(heatmap_data)[1:16],
             col = c("blue1", "skyblue", "white", "pink", "red"))


module.gene.mapping <- as.data.frame(all_subtype_bwnet$colors)
module <- rownames(module.gene.mapping)[module.gene.mapping$`all_subtype_bwnet$colors` == "turquoise"]

temp <- approved_openTargets[approved_openTargets$`Target ID` %in% module, ]
table(unique(approved_openTargets$`Target Approved Symbol`) %in% temp$`Target Approved Symbol`)


gene.signf.corr <- cor(all_wgcna_data, traits$data.tumour.vs.all, use = 'p')
gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples = nrow(all_wgcna_data))
gene.signf.corr.pvals <- rownames_to_column(as.data.frame(gene.signf.corr.pvals))
gene.signf.corr.pvals <- gene.signf.corr.pvals[order(-gene.signf.corr.pvals$V1), ]
rownames(gene.signf.corr.pvals) <- NULL

temp <- gene.signf.corr.pvals[gene.signf.corr.pvals$rowname %in% unique(approved_openTargets$`Target ID`), ]
temp <- merge(temp, approved_openTargets, by.x = "rowname", by.y = "Target ID")






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










