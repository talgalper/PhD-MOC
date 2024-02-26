library(biomaRt)
library(ggplot2)
library(edgeR)
library(reshape2)
library(ggbreak)
options(scipen = 999)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene_id <- getBM(attributes = c("external_gene_name", "ensembl_gene_id"), 
                 filters = "external_gene_name", 
                 values = c("ESR2", "ESR1", "PGR", "ERBB2", "MKI67"), 
                 mart = ensembl)



load("RData/lumA/BRCA_lumA_DE_data.RData")
load("RData/lumB/BRCA_lumB_DE_data.RData")
load("RData/Her2/BRCA_Her2_DE_data.RData")
load("RData/Basal/BRCA_Basal_DE_data.RData")
load("RData/Normal/BRCA_Normal_DE_data.RData")

# logCPM on raw data
LumA_CPM <- as.data.frame(cpm(LumA_unstranded))
LumB_CPM <- as.data.frame(cpm(LumB_unstranded))
Her2_CPM <- as.data.frame(cpm(Her2_unstranded))
Basal_CPM <- as.data.frame(cpm(Basal_unstranded))
TCGA_normal_CPM <- as.data.frame(cpm(TCGA_normal_unstranded))
GTEx_normal_CPM <- as.data.frame(cpm(GTEx_data))

# subset target genes
LumA_subset <- LumA_CPM[rownames(LumA_CPM) %in% gene_id$ensembl_gene_id, ]
LumB_subset <- LumB_CPM[rownames(LumB_CPM) %in% gene_id$ensembl_gene_id, ]
Her2_subset <- Her2_CPM[rownames(Her2_CPM) %in% gene_id$ensembl_gene_id, ]
Basal_subset <- Basal_CPM[rownames(Basal_CPM) %in% gene_id$ensembl_gene_id, ]
GTEx_normal_subset <- GTEx_normal_CPM[rownames(GTEx_normal_CPM) %in% gene_id$ensembl_gene_id, ]
TCGA_normal_subset <- TCGA_normal_CPM[rownames(TCGA_normal_CPM) %in% gene_id$ensembl_gene_id, ]

# create function for SE
std.error <- function(x) sd(x)/sqrt(length(x))

# apply function over CPM data
SE_df <- data.frame(LumA = apply(LumA_subset, 1, std.error),
                    LumB = apply(LumB_subset, 1, std.error),
                    Her2 = apply(Her2_subset, 1, std.error),
                    Basal = apply(Basal_subset, 1, std.error),
                    normal_GTEx = apply(GTEx_normal_subset, 1, std.error),
                    normal_TCGA = apply(TCGA_normal_subset, 1, std.error))

# calculate mean of CPM for each gene
mean_CPM_df <- data.frame(LumA = apply(LumA_subset, 1, mean),
                            LumB = apply(LumB_subset, 1, mean),
                            Her2 = apply(Her2_subset, 1, mean),
                            Basal = apply(Basal_subset, 1, mean),
                            normal_GTEx = apply(GTEx_normal_subset, 1, mean),
                            normal_TCGA = apply(TCGA_normal_subset, 1, mean))

# put gene IDs back in rownames
gene_id <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), 
                                     filters = "ensembl_gene_id", 
                                     values = rownames(mean_CPM_df), 
                                     mart = ensembl)
rownames(mean_CPM_df) <- gene_id$external_gene_name
rownames(SE_df) <- gene_id$external_gene_name


# Transpose the data for plotting
mean_CPM_df <- t(mean_CPM_df)
mean_CPM_df <- as.data.frame(mean_CPM_df)

SE_df <- t(SE_df)
SE_df <- as.data.frame(SE_df)

# Add gene names as a column
mean_CPM_df$subtype <- rownames(mean_CPM_df)
SE_df$subtype <- rownames(SE_df)

# Melt the data for plotting
mean_CPM_df <- melt(mean_CPM_df, id.vars = "subtype")
colnames(mean_CPM_df)[2] <- "gene"
colnames(mean_CPM_df)[3] <- "mean_cpm"


SE_df <- melt(SE_df, id.vars = "subtype")
colnames(SE_df)[2] <- "gene"
colnames(SE_df)[3] <- "SE"

# combine data
plot_data <- merge(mean_CPM_df, SE_df)


ggplot(plot_data, aes(x = gene, y = mean_cpm, fill = subtype)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = mean_cpm - SE, ymax = mean_cpm + SE), position = position_dodge(width = 0.8), width = 0.25) +
  labs(x = "Gene", y = "Mean Expression Level (CPM)", fill = "Cancer Subtype") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_fill_brewer(palette = "Paired") +
  scale_y_break(c(1000, 4200))










