library(data.table)
library(tidyverse)
library(WGCNA)
library(edgeR)
library(DESeq2)

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


sft <- pickSoftThreshold(all_wgcna_data,
                         blockSize = 45000,
                         verbose = 2)

sft <- sft$fitIndices

library(gridExtra)
library(grid)
a1 <- ggplot(sft, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1, size = 7) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit') +
  theme_classic() +
  theme(axis.text = element_text(size = 18, colour = "black"),
        axis.title = element_text(size = 20))

a2 <- ggplot(sft, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 1000, size = 7) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic() +
  theme(axis.text = element_text(size = 18, colour = "black"),
        axis.title = element_text(size = 20))

grid.arrange(a1, a2, nrow = 2)
rm(a1, a2)
# separate and create dissTOM's
tumour_data <- all_subtype_counts_filt[, colnames(all_subtype_counts_filt) %in% sample_info$sample[sample_info$group != "control"]]
control_data <- all_subtype_counts_filt[, colnames(all_subtype_counts_filt) %in% sample_info$sample[sample_info$group == "control"]]

tumour_TOM <- adjacency(tumour_data, power = 12, type = "signed")

## perform diff_i method

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




