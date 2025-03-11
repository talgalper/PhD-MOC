library(data.table)
library(tidyverse)
library(WGCNA)
library(edgeR)
library(DESeq2)
library(doParallel)

registerDoParallel(cores = 30)
enableWGCNAThreads(nThreads = 30)
WGCNAnThreads()

## create WGCNA adjacencey matrices
# load in data
load("../BRCA_pipe/RData/LumA/DE_data.RData")
load("../BRCA_pipe/RData/LumB/DE_data.RData")
load("../BRCA_pipe/RData/Her2/DE_data.RData")
load("../BRCA_pipe/RData/basal/DE_data.RData")
load("../BRCA_pipe/RData/TCGA_normal.RData")

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

# add sample type to sample_info
load("../BRCA_pipe/RData/TCGA_query.RData")
common <- common[, c(1,3)]
sample_info <- merge(sample_info, common, by.x = "sample", by.y = "cases", all.x = T)
sample_info$sample_type <- ifelse(is.na(sample_info$sample_type), "Healthy", sample_info$sample_type)

# combine all tumour samples
all_subtypes <- cbind(LumA_unstranded, LumB_unstranded, Her2_unstranded, Basal_unstranded)

# clean env
rm(normal_unstranded, LumA_unstranded, LumB_unstranded, Her2_unstranded, Basal_unstranded)
rm(control_info, lumA_info, lumB_info, her2_info, basal_info)
rm(clinical, common, query_TCGA, subtypes)
collectGarbage()

STN_samples <- sample_info$sample[sample_info$sample_type == "Solid Tissue Normal"]
all_subtypes <- all_subtypes[, !colnames(all_subtypes) %in% STN_samples]

# QC + combines tumour and control samples
all_subtype_counts_filt <- filter_low_expr(tumour_matrix = all_subtypes,
                                           control_matrix = GTEx_ENS)

# normalisation (transposes matrix)
all_wgcna_data <- vst_norm(all_subtype_counts_filt)


# separate and create dissTOM's
tumour_data <- all_subtype_counts_filt[, colnames(all_subtype_counts_filt) %in% sample_info$sample[sample_info$group != "control"]]
control_data <- all_subtype_counts_filt[, colnames(all_subtype_counts_filt) %in% sample_info$sample[sample_info$group == "control"]]

tumour_data[] <- lapply(tumour_data, as.numeric)
tumour_data <- t(tumour_data)
tumour_TOM <- TOMsimilarityFromExpr(tumour_data, TOMType = "unsigned", power = 6, nThreads = 30)

control_data[] <- lapply(control_data, as.numeric)
control_data <- t(control_data)
control_TOM <- TOMsimilarityFromExpr(control_data, TOMType = "unsigned", power = 6, nThreads = 30)

## perform diff_i method
tumour_dissTOM <- 1 - tumour_TOM


## figure out how they performed wierd STRING PPI corss reference

## create index and indexed edge list


# create score file
load("../BRCA_pipe/latest_run/RData/STN_filt/dif_exp.RData")
DE_data <- subset(dif_exp, select = c("gene_id", "logFC"))
DE_data$logFC_abs <- abs(DE_data$logFC) # get absolute values
DE_data_abs <- subset(DE_data, select = c("gene_id", "logFC_abs"))


fwrite(gene_index, "BRCA/STN_filt/data/gene_index.tsv", col.names = F, sep = "\t")
fwrite(edge_list_index, "BRCA/STN_filt/data/edge_list_index.tsv", col.names = F, sep = "\t")
fwrite(DE_data_abs, "BRCA/STN_filt/data/updated_DE/logFC_scores_abs.tsv", col.names = F, sep = "\t")




