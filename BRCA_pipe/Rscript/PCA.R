library(tidyverse)
library(edgeR)

# Load data
load("RData/LumA/DE_data.RData")
load("RData/LumB/DE_data.RData")
load("RData/Her2/DE_data.RData")
load("RData/basal/DE_data.RData")
load("RData/TCGA_normal.RData")

GTEx_data <- read.table("gene_reads_2017-06-05_v8_breast_mammary_tissue.gct", skip = 2, header = T)
GTEx_data <- GTEx_data[, -1]
rownames(GTEx_data) <- NULL
GTEx_data$Name <- gsub("\\.\\d+", "", GTEx_data$Name)

# opt for having gene Ensembl IDs instead of gene names as rownames (same as TCGA)
GTEx_ENS <- column_to_rownames(GTEx_data, "Name")
GTEx_ENS <- GTEx_ENS[ , -1]
rownames <- rownames(GTEx_ENS)
GTEx_ENS <- as.data.frame(sapply(GTEx_ENS, as.numeric))
rownames(GTEx_ENS) <- rownames
rm(GTEx_data, rownames)


## all data in one PCA
TCGA_data <- cbind(normal_unstranded, LumA_unstranded, LumB_unstranded, Her2_unstranded, Basal_unstranded)
data_master <- merge(GTEx_ENS, TCGA_data, by = "row.names")
data_master <- column_to_rownames(data_master, "Row.names")

sample_info <- data.frame(
  sample = colnames(data_master),
  group = c(rep("GTEx", ncol(GTEx_ENS)), 
            rep("TCGA_normal", ncol(normal_unstranded)),
            rep("LumA", ncol(LumA_unstranded)),
            rep("LumB", ncol(LumB_unstranded)),
            rep("Her2", ncol(Her2_unstranded)),
            rep("Basal", ncol(Basal_unstranded))))
sample_info$group <- factor(sample_info$group, levels = c("GTEx", "TCGA_normal", "LumA", "LumB", "Her2", "Basal"))

data_master_filt <- filterByExpr(data_master, group = sample_info$group)
data_master_filt <- data_master[data_master_filt, ]

data_master_norm <- cpm(as.matrix(data_master_filt), log = T)

# requires function from WGCNA_functions.R
PCA_all_data <- plot_PCA(expr_data = t(data_master_filt), sample_info = sample_info)




## tumour and control data separately
tumour_data <- cbind(LumA_unstranded, LumB_unstranded, Her2_unstranded, Basal_unstranded)

sample_info <- data.frame(
  sample = colnames(tumour_data),
  group = c(rep("LumA", ncol(LumA_unstranded)),
            rep("LumB", ncol(LumB_unstranded)),
            rep("Her2", ncol(Her2_unstranded)),
            rep("Basal", ncol(Basal_unstranded))))
sample_info$group <- factor(sample_info$group, levels = c("LumA", "LumB", "Her2", "Basal"))

tumour_filt <- filterByExpr(tumour_data, group = sample_info$group)
tumour_filt <- tumour_data[tumour_filt, ]

tumour_norm <- cpm(as.matrix(tumour_filt), log = T)

tumour_PCA <- plot_PCA(expr_data = t(logCPM), sample_info = sample_info, plot_tree = F)





control_data <- merge(GTEx_ENS, normal_unstranded, by = "row.names")
control_data <- column_to_rownames(control_data, "Row.names")

sample_info <- data.frame(
  sample = colnames(control_data),
  group = c(rep("GTEx", ncol(GTEx_ENS)),
            rep("Normal_unstranded", ncol(normal_unstranded))))
sample_info$group <- factor(sample_info$group, levels = c("GTEx", "Normal_unstranded"))

control_filt <- filterByExpr(control_data, group = sample_info$group)
control_filt <- control_data[control_filt, ]

control_norm <- cpm(as.matrix(control_filt), log = T)
control_PCA <- plot_PCA(expr_data = t(control_norm), sample_info = sample_info)





tumour_data <- list(LumA_unstranded = LumA_unstranded, 
                    LumB_unstranded = LumB_unstranded, 
                    Her2_unstranded = Her2_unstranded, 
                    Basal_unstranded = Basal_unstranded)

# form filtering for each group separately
tumour_data_filt <- list() 
for (i in seq_along(tumour_data)) {
  subtype <- names(tumour_data)[i]
  data <- tumour_data[[i]]
  cat("processing", subtype, "\n")
  
  filt <- filterByExpr(data)
  filt <- data[filt, ]
  
  tumour_data_filt[[subtype]] <- filt
  
  rm(subtype, data, filt, i)
}


tumour_filt <- filterByExpr(tumour_data)
tumour_filt <- merged_df[tumour_filt, ]



# combine by common genes
tumour_data <- cbind(LumA_unstranded, LumB_unstranded, Her2_unstranded, Basal_unstranded)
common_genes <- rownames(tumour_data)[rownames(tumour_data) %in% rownames(GTEx_ENS)]
tumour_data <- tumour_data[rownames(tumour_data) %in% common_genes, ]
