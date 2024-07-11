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

hist(as.matrix(tumour_norm_counts))
hist(as.matrix(tumour_exp_counts))


# QC + combines tumour and control samples
all_subtype_counts_filt <- filter_low_expr(tumour_matrix = tumour_exp_counts,
                                           control_matrix = GTEx_exp_counts,
                                           sep = T)

hist(as.matrix(all_subtype_counts_filt$tumour))


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