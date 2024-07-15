library(tidyverse)
library(WGCNA)

norm_counts <- read.table("../../../../Downloads/TCGA-GTEx-TARGET-gene-exp-counts.deseq2-normalized.log2", header = T, row.names = 1)
exp_counts <- read.table("../../../../Downloads/TcgaTargetGtex_gene_expected_count", header = T, row.names = 1)

phenotype_data <- read.table("../../../../Downloads/TcgaTargetGTEX_phenotype.txt", sep = "\t", header = T)
phenotype_data$sample <- gsub("-", ".", phenotype_data$sample)

BRCA_sample_info <- phenotype_data[phenotype_data$X_primary_site == "Breast", ]

tumour_info <- BRCA_sample_info[BRCA_sample_info$X_sample_type == "Primary Tumor", ]
GTEx_info <- BRCA_sample_info[BRCA_sample_info$X_sample_type == "Normal Tissue", ]
TCGA_normal_info <- BRCA_sample_info[BRCA_sample_info$X_sample_type == "Solid Tissue Normal", ]

tumour_exp_counts <- exp_counts[, colnames(exp_counts) %in% tumour_info$sample]
GTEx_exp_counts <- exp_counts[, colnames(exp_counts) %in% GTEx_info$sample]
TCGA_normal_exp_counts <- exp_counts[, colnames(exp_counts) %in% TCGA_normal_info$sample]

tumour_norm_counts <- norm_counts[, colnames(norm_counts) %in% tumour_info$sample]
GTEx_norm_counts <- norm_counts[, colnames(norm_counts) %in% GTEx_info$sample]
TCGA_normal_norm_counts <- norm_counts[, colnames(norm_counts) %in% TCGA_normal_info$sample]

save(tumour_info, GTEx_info, TCGA_normal_info, phenotype_data,
     tumour_exp_counts, GTEx_exp_counts, TCGA_normal_exp_counts,
     tumour_norm_counts, GTEx_norm_counts, TCGA_normal_norm_counts,
     file = "../../../../Desktop/xena_data.RData")


load("../../../../Desktop/xena_data.RData")

hist(as.matrix(tumour_norm_counts))
hist(as.matrix(tumour_exp_counts))


# QC + combines tumour and control samples
library(edgeR)
all_subtype_counts_filt <- filter_low_expr(tumour_matrix = tumour_norm_counts,
                                           control_matrix = GTEx_norm_counts,
                                           sep = T)

hist(as.matrix(all_subtype_counts_filt$tumour))
hist(as.matrix(all_subtype_counts_filt$control))


tumour_wgcna_data <- t(all_subtype_counts_filt$tumour)
control_wgcna_data <- t(all_subtype_counts_filt$control)

# normalisation (transposes matrix)
#tumour_wgcna_data <- vst_norm(all_subtype_counts_filt$tumour)

# plot PCA
control_info <- data.frame(sample = colnames(all_subtype_counts_filt$control),
                           group = rep("control", ncol(all_subtype_counts_filt$control)))
tumour_info <- data.frame(sample = colnames(all_subtype_counts_filt$tumour),
                        group = rep("tumour", ncol(all_subtype_counts_filt$tumour)))

PCA_results_tumour <- plot_PCA(expr_data = tumour_wgcna_data,
                             sample_info = tumour_info,
                             plot_tree = F,
                             output_plot_data = T)

PCA_results_control <- plot_PCA(expr_data = control_wgcna_data,
                             sample_info = control_info,
                             plot_tree = F,
                             output_plot_data = T)



# choose soft thresholding power
tumour_sft_data <- pick_power(WGCNA_data = tumour_wgcna_data,
                              network_type = "unsigned")

control_sft_data <- pick_power(WGCNA_data = control_wgcna_data,
                               network_type = "unsigned")



# RESTART R AND LOAD WGCNA ONLY
library(WGCNA)
library(doParallel)
nCores = 8
registerDoParallel(cores = nCores)
enableWGCNAThreads(nThreads = nCores)
WGCNAnThreads()

control_bwnet <- network_modules(WGCNA_data = control_wgcna_data,
                                 Power = 8)
tumour_bwnet <- network_modules(WGCNA_data = tumour_wgcna_data,
                                Power = 8)

save(control_bwnet, file = "BRCA/RData/GTEx/control_bwnet.RData")
save(tumour_bwnet, file = "BRCA/RData/GTEx/tumour_bwnet.RData")






# module preservation using expr data
common_genes <- intersect(colnames(tumour_wgcna_data), colnames(control_wgcna_data))
tumour_common <- tumour_wgcna_data[, colnames(tumour_wgcna_data) %in% common_genes]
control_common <- control_wgcna_data[, colnames(control_wgcna_data) %in% common_genes]

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



plot_data <- modulePreservation_plt$plot_data$plot_data
non_preserved_modules <- plot_data[plot_data$medianRank.pres < 8 & plot_data$Zsummary.pres > 10, ]

nonPreservedGenes <- common_genes[tumour_bwnet$colors %in% non_preserved_modules$cluster]
nonPreservedGenes <- gsub("\\.\\d+", "", nonPreservedGenes)


# read in drug data
OpenTargets <- read.csv("../BRCA_pipe/OpenTargets_data/OpenTargets_unique_drug.csv", row.names = 1)
OpenTargets_raw <- read_tsv("../BRCA_pipe/OpenTargets_data/breast_carcinoma_known_drugs.tsv")

NIH_targets <- read.table("../BRCA_pipe/NIH_BRCA_approved_drugs.txt", sep = "\t")
colnames(NIH_targets)[1] <- "approved_drugs"
NIH_targets$approved_drugs <- toupper(NIH_targets$approved_drugs)

approved_openTargets <- merge(NIH_targets, OpenTargets_raw, by.x = "approved_drugs", by.y = "Drug Name")
approved_openTargets <- approved_openTargets[!duplicated(approved_openTargets$approved_drugs) | !duplicated(approved_openTargets$`Target ID`), ]


temp <- approved_openTargets[approved_openTargets$`Target ID` %in% nonPreservedGenes, ]

table(unique(approved_openTargets$`Target Approved Symbol`) %in% unique(temp$`Target Approved Symbol`))
table(unique(OpenTargets$Target.ID) %in% nonPreservedGenes)



colours <- labels2colors(tumour_bwnet$colors)
tumour_kWithin <- intramodularConnectivity.fromExpr(tumour_data, colours, power = 6)
rownames(tumour_kWithin) <- colnames(tumour_data)








