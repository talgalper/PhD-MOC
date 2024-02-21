library(biomaRt)
library(ggplot2)
library(edgeR)
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

LumA_expression <- LumA_unstranded[rownames(LumA_unstranded) %in% gene_id$ensembl_gene_id, ]
LumB_expression <- LumB_unstranded[rownames(LumB_unstranded) %in% gene_id$ensembl_gene_id, ]
Her2_expression <- Her2_unstranded[rownames(Her2_unstranded) %in% gene_id$ensembl_gene_id, ]
Basal_expression <- Basal_unstranded[rownames(Basal_unstranded) %in% gene_id$ensembl_gene_id, ]
GTEx_normal_expression <- GTEx_data[rownames(GTEx_data) %in% gene_id$ensembl_gene_id, ]
TCGA_normal_expression <- TCGA_normal_unstranded[rownames(TCGA_normal_unstranded) %in% gene_id$ensembl_gene_id, ]

expression_df <- data.frame(LumA = apply(LumA_expression, 1, mean),
                            LumB = apply(LumB_expression, 1, mean),
                            Her2 = apply(Her2_expression, 1, mean),
                            Basal = apply(Basal_expression, 1, mean),
                            normal_GTEx = apply(GTEx_normal_expression, 1, mean),
                            normal_TCGA = apply(TCGA_normal_expression, 1, mean))

gene_id <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), 
                                     filters = "ensembl_gene_id", 
                                     values = rownames(expression_df), 
                                     mart = ensembl)

rownames(expression_df) <- gene_id$external_gene_name


# Transpose the data for plotting
expression_df <- t(expression_df)
expression_df <- as.data.frame(expression_df)

# Add gene names as a column
expression_df$Gene <- rownames(expression_df)

# Melt the data for plotting
expression_df <- melt(expression_df, id.vars = "Gene")

cpm_expression <- expression_df
cpm_expression$value <- cpm(expression_df$value, log = T)
colnames(cpm_expression)[3] <- "value"

ggplot(expression_df, aes(x = variable, y = value, fill = Gene)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  labs(x = "Gene", y = "Expression Level (counts)", fill = "Cancer Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_fill_brewer(palette = "Paired")

ggplot(cpm_expression, aes(x = variable, y = value, fill = Gene)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  labs(x = "Gene", y = "Expression Level (logCPM)", fill = "Cancer Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_fill_brewer(palette = "Paired")



